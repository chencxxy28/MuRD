#create the test data set
testset<-function(eset,list.marker)
{
    bulk.eset<-eset
    #bulk_data<-(exprs(bulk.eset))
    marker_gene_used<-unlist(list.marker)
    gene_test<-setdiff(rownames(bulk.eset),marker_gene_used)
    bulk_data<-bulk.eset[gene_test,]
    bulk_data
}

#A criterion function to assess the performance of deconvolution
criteria.onegroup<-function (bulk.data, prop.used)
{
    message("calculate criteria")
    bulk_matrix_raw<-(exprs(bulk.data))
    find_zero_index<-which(rowSums(bulk_matrix_raw,na.rm=T)*1e5==0)
    if(length(find_zero_index)>1)
    {
        bulk_nozero<-getCPM0(bulk_matrix_raw[-find_zero_index,])
    }else{
        bulk_nozero<-getCPM0(bulk_matrix_raw)
    }
    all_nozero_index<-which(apply(bulk_nozero,1,function (x) sum(x==0))==0)
    bulk_all_nozero<-bulk_nozero[all_nozero_index,]
    check_rowsum<-rowSums(bulk_all_nozero)
    bulk_all_nozero<-bulk_all_nozero[check_rowsum<quantile(check_rowsum,0.95) & check_rowsum>quantile(check_rowsum,0.15),]
    sigma_subject<-apply(bulk_all_nozero,1,sd)
    bulk_all_nozero<-bulk_all_nozero[sigma_subject!=0,]
    gene_validate<-rownames(bulk_all_nozero)


    prop_new<-as.matrix(prop.used)
    ##do training

    bulk_nozero_train<-bulk_all_nozero
    prop_new_train<-prop_new
    X_all<-t(sapply(1:length(gene_validate),function (xx)
    {
        gene_xx<-gene_validate[xx]
        bulk_xx<-bulk_nozero_train[gene_xx,]*1e5
        prop_xx<-prop_new_train
        fit<-nnls(A=prop_xx,B=bulk_xx)
        fit$X
    }
    ))
    rownames(X_all)<-gene_validate

    bulk_nozero_test<-bulk_all_nozero
    prop_new_test<-prop_new
    prop_new_withscaler<-t(sapply(1:ncol(bulk_nozero_test), function (x)
    {
        #each subject estimation
        Y_initial<-bulk_nozero_test[,x]
        zero_index<-which(Y_initial==0) #find genes with zero expression in current subject
        if(length(zero_index)>0)
        {
            Y_initial<-Y_initial[-zero_index] #delate zero expressed genes
        }
        genes_left<-names(Y_initial)

        #construct Y
        gene_used_now<-intersect(genes_left,gene_validate)
        Y<-Y_initial[gene_used_now]*1e5
        #Y<-Y[Y<quantile(Y,0.85) & Y>quantile(Y,0.15)] #remove the outliers
        gene_used_final<-names(Y)

        #construct X
        X<-X_all[gene_used_final,]

        values_final<-sort(abs(Y-X%*%(prop_new_test[x,])),decreasing = T)[1:50]

        c(median(values_final),mean(values_final),sqrt(mean(values_final^2)))
        #print(x)
    })
    )
    apply(prop_new_withscaler,2,mean)
}

#extract estimated cell type proportions via MuRD
MuRD.predict.prop<-function(murd.output)
{
    goodnessoffit<-murd.output$metrics
    murd.output$est[[which.min(goodnessoffit)]]
}

#extract estimated cell type proportions via ref-free
reffree.predict.prop<-function(murd.output)
{
    murd.output$est[[1]]
}

#evaluation
evaluate<-function(est.prop,true.prop)
{
    ct.mad<-apply(est.prop-true.prop,2,function (x) mean(abs(x)))
    all.mad<-mean(ct.mad)
    ct.ccc<-unlist(sapply(1:ncol(est.prop),function(x)
    {
        CCC_all<-CCC(est.prop[,x],true.prop[,x])
        CCC_all$rho[1]
    }
        )
    )
    all.ccc<-mean(ct.ccc)
    ct.cor<-unlist(sapply(1:ncol(est.prop),function(x)
    {
        cor(est.prop[,x],true.prop[,x])
    }
    )
    )
    all.cor<-mean(ct.cor)
    cell.type.eva<-data.frame(ct.mad=ct.mad,ct.ccc=ct.ccc,ct.cor=ct.cor)
    all.eva<-data.frame(all.mad=all.mad, all.ccc=all.ccc,all.cor=all.cor)
    return(list(cell.type.eva=cell.type.eva,all.eva=all.eva))
}

#to obtain population-level cell type proportions from a given subject-level cell type proportions
pop.ct.prop.subj<-function(subj.prop)
{
    apply(subj.prop,2,mean)
}

#to obtain population-level cell type proportions from single cell RNA-seq data
pop.ct.prop.scRNA<-function(scRNA,cluster,sample,sep.group=F, group=NULL)
{
    ct_ad<-unique(scRNA@phenoData@data[[cluster]])
    sampleid<-unique(scRNA@phenoData@data[[sample]])

    sc_proportions<-sapply(1:length(unique(scRNA@phenoData@data[[sample]])),function (x)
    {
        ct_x<-scRNA@phenoData@data[[cluster]][sampleid==sampleid[x]]
        table(ct_x)/sum(table(ct_x))
    }
    )
    sc_proportions<-t(sc_proportions[ct_ad,])
    rownames(sc_proportions)<-sampleid
    if (sep.group == T)
    {
        group_subid<-tapply(scRNA@phenoData@data[[group]], scRNA@phenoData@data[[sample]],unique)
        group_label<-unique(group_subid)
        output<-t(sapply(1:(length(group_label)),function (x){
             apply(sc_proportions[names(which(group_label==group_label[x])),],2,mean)
         }))
        colnames(output)<-ct_ad
        rownames(output)<-group_label
    }else{
        output<-apply(sc_proportions,2,mean)
    }
    return(pop.ct.prop=output)
}

