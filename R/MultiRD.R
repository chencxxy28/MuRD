#estimate mean profile for each marker gene, given cell proportions
estimate.geneprofile<-function(bulk.data,gene.used,celltype.unique,cluster.identifier,prop.est)
{
    X<-t(sapply(1:length(gene.used),function (xx)
    {
        x_xx<-rep(0,length(celltype.unique))
        gene_xx<-gene.used[xx]
        bulk_xx<-bulk.data[gene_xx,]
        cell_type_index<-cluster.identifier[xx]
        prop_xx<-prop.est[,cell_type_index]
        #x_xx[cell_type_index]<-sum((bulk_xx*prop_xx))/sum(prop_xx^2)*1e5
        x_xx[cell_type_index]<-sum(bulk_xx)/sum(prop_xx)*1e5 #this formula is weird, I think the previous one is better
        x_xx
    }
    ))
    X
}

#iterating updating cell type proportions and gene mean profiles
MultiRD.onegroup<-function (bulk.data,list.marker,celltype.unique,subject.level.proportion,population.level.proportion,proportion.sd=1,lambda.option=c(seq(from=0,to=0.075,length=15),10,50,100,500,1000),tol=0.001,iter.num=1000){
    #construct test set for evaluation
    #bulk.eset<-bulk.data
    #bulk_data<-(exprs(bulk.eset))
    #marker_gene_used<-unlist(list.marker)
    #gene_test<-setdiff(rownames(bulk_data),marker_gene_used)
    #bulk_data<-bulk_data[gene_test,]
    #exprs_data<-as.matrix(bulk_data)
    #pdata<-data.frame(sample=colnames(bulk_data))
    #fdata<-data.frame(genes=rownames(bulk_data))
    #rownames(pdata)<-colnames(bulk_data)
    #rownames(fdata)<-rownames(bulk_data)
    #bulk.data_test<-ExpressionSet(exprs_data,
    #                                 AnnotatedDataFrame(pdata),
    #                                 AnnotatedDataFrame(fdata))
    bulk.data_test<-testset(eset=bulk.data,list.marker=list.marker)

    #run reference free approach
    ct.sub<-celltype.unique
    bulk.eset<-bulk.data
    bulk_matrix_raw<-(exprs(bulk.eset))
    find_zero_index<-which(rowSums(bulk_matrix_raw,na.rm=T)*1e5==0)
    if(length(find_zero_index)>0)
    {
        bulk_nozero<-getCPM0(bulk_matrix_raw[-find_zero_index,])
    }else{
        bulk_nozero<-getCPM0(bulk_matrix_raw)
    }
    marker_all<-unlist(list.marker)
    #lambda.option<-c(seq(from=0,to=0.075,length=15),1,10,30,50,1000)
    #lambda.option<-c(1)
    #group<-bulk.data@phenoData@data[["groups"]]

    marker_list_sub<-lapply(1:length(list.marker),function (xx)
    {
        marker_xx<-intersect(list.marker[[xx]],rownames(bulk_nozero)) #create new marker list
    }
    )
    names(marker_list_sub)<-ct.sub
    gene_used<-unlist(marker_list_sub) #vectorize the values
    cluster_identifier<-unlist(lapply(1:length(marker_list_sub),function (xx)
    {
        rep(xx,length(marker_list_sub[[xx]])) #create the identifier to discriminate which cell type are markers coming from
    }
    ))

    criterion<-NULL
    est_all<-list(NULL)

    for(ll in 1:length(lambda.option))
    {
        iter<-0
        prop_old<-matrix(1/length(ct.sub),nrow=ncol(bulk_matrix_raw),ncol=length(ct.sub))
        lambda<-lambda.option[ll]
        repeat{
            #update mean gene profile
            X<-estimate.geneprofile(bulk.data=bulk_nozero,gene.used=gene_used,celltype.unique=ct.sub,cluster.identifier=cluster_identifier,prop.est=prop_old)
            rownames(X)<-gene_used

            prop_new<-t(sapply(1:ncol(bulk_nozero), function (x)
            {
                #each subject estimation
                Y_initial<-bulk_nozero[,x]
                zero_index<-which(Y_initial==0) #find genes with zero expression in current subject
                if (length(zero_index)>0)
                {
                    Y_initial<-Y_initial[-zero_index] #delate zero expressed genes !!!!!!!!!!!!check this!!!!!!!!!!!
                }
                #Y_initial<-Y_initial[Y_initial<quantile(Y_initial,0.99)] #remove the outliers
                gene_names<-names(Y_initial)
                marker_genes<-intersect(gene_names,marker_all) #match the marker gene
                #check marker gene availability for each cell types
                marker_list_sub_i<-lapply(1:length(list.marker),function (xx)
                {
                    marker_xx<-intersect(list.marker[[xx]],marker_genes) #create new marker list
                }
                )
                names(marker_list_sub_i)<-ct.sub
                gene_used<-unlist(marker_list_sub_i) #vectorize the values
                cluster_identifier<-unlist(lapply(1:length(marker_list_sub_i),function (xx)
                {
                    rep(xx,length(marker_list_sub_i[[xx]])) #create the identifier to discriminate which cell type are markers coming from
                }
                ))
                find_zero_counts<-which(table(cluster_identifier)==0)

                #construct Y
                Y<-Y_initial[gene_used]*1e5
                lambda_adjust<-sum(max(Y)^2)*lambda
                if(ll>=length(lambda.option)-4)
                {
                    ave_truep<-as.matrix(subject.level.proportion)[x,]
                }else{
                    ave_truep<-population.level.proportion
                }

                Y_aug<-ave_truep*sqrt(lambda_adjust)/proportion.sd
                #Y_aug<-c(sapply(1:length(combo), function (xx) {sqrt(lambda_adjust*weight[xx])*(combo[[xx]][x,])}))
                Y_comb<-c(Y,Y_aug)

                #construct X matrix
                #X<-t(sapply(1:length(gene_used),function (xx)
                #{
                #     x_xx<-rep(0,length(ct.sub))
                #     gene_xx<-gene_used[xx]
                #     bulk_xx<-bulk_nozero[gene_xx,]
                #     cell_type_index<-cluster_identifier[xx]
                #     prop_xx<-prop_old[,cell_type_index]
                #     #x_xx[cell_type_index]<-sum((bulk_xx*prop_xx))/sum(prop_xx^2)*1e5
                #     x_xx[cell_type_index]<-sum(bulk_xx)/sum(prop_xx)*1e5 #this formula is weird, I think the previous one is better
                #     x_xx
                # }
                # ))
                X<-X[gene_used,]
                X_aug<-sqrt(lambda_adjust)*diag(1,length(ct.sub))/proportion.sd

                #X_aug<-do.call(rbind,lapply(1:length(combo), function (xx) {sqrt(lambda_adjust*weight[xx])*diag(1,length(ct.sub))}))
                X_comb<-rbind(X,X_aug)
                #heatmap(X)
                ##construct the constrain matrix
                E_used<-rep(1,length(ct.sub))
                F_used<-1
                G_used<-diag(1,length(ct.sub))
                H_used<-rep(0,length(ct.sub))
                #fit<-lsei(A=X,B=Y,E=E_used,F=F_used,G=G_used,H=H_used)
                #fit$X
                fit<-nnls(A=X_comb,B=Y_comb)

                #calcualte GCV values
                #x_x<-t(X)%*%X
                #inverse_term<-solve(x_x+sum(lambda_adjust*weight))
                #P_lambda<-X%*%inverse_term%*%t(X)
                #penalty_1<-(1-(tr(P_lambda)-1)/length(Y))^2
                #penalty_2<-
                #gcv_values<-(Y%*%((diag(1,nrow(P_lambda))-P_lambda)^2)%*%Y)/penalty_1+
                fit$X/sum(fit$X)
                #print(x)
            })
            )
            if(mean(abs(prop_new[,1:length(ct.sub)]-prop_old))<tol | iter>iter.num)
            {
                break
            }else{
                prop_old<-prop_new[,1:length(ct.sub)]
                iter<-iter+1
            }
        }
        est_meta_all<-prop_new
        criterion_sub<-criteria.onegroup(bulk.data=bulk.data_test, prop.used=prop_new)[2]
        criterion<-c(criterion,criterion_sub)
        est_all[[ll]]<-est_meta_all
        print(ll)
    }

    prop_est<-est_all[[which.min(criterion)]]
    list(est=est_all,metrics=criterion)
}
