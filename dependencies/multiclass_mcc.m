function mcc=multiclass_mcc(cm_svm_array)
%INPUT- cm_svm_array is a confusion matrix
%OUTPUT- MCC for multiclass confusion matrix
% https://stats.stackexchange.com/questions/187768/matthews-correlation-coefficient-with-multi-class
mcc_numerator=0;count=1;
% limits klm=1 TO n SUM(ckk.cml - clk.ckm)
for k = 1:1:length(cm_svm_array)
    for l=1:1:length(cm_svm_array)
        for m=1:1:length(cm_svm_array)
          mcc_numerator1(count) = (cm_svm_array(k,k) *cm_svm_array(m,l))-...
                                  (cm_svm_array(l,k)*cm_svm_array(k,m));
          mcc_numerator=mcc_numerator+mcc_numerator1(count);
          count=count+1;
        end
    end
end

mcc_denominator_1=0 ; count=1;
for k=1:1:length(cm_svm_array)
     mcc_den_1_part1=0;
    for l=1:1:length(cm_svm_array)
        mcc_den_1_part1= mcc_den_1_part1+cm_svm_array(l,k);
    end
    mcc_den_1_part2=0;
    for f=1:1:length(cm_svm_array)
        if f ~=k
          for g=1:1:length(cm_svm_array)
            mcc_den_1_part2= mcc_den_1_part2+cm_svm_array(g,f);
          end
        end
    end
    mcc_denominator_1=(mcc_denominator_1+(mcc_den_1_part1*mcc_den_1_part2));
end

mcc_denominator_2=0; count=1;
for k=1:1:length(cm_svm_array)
     mcc_den_2_part1=0;
    for l=1:1:length(cm_svm_array)
        mcc_den_2_part1= mcc_den_2_part1+cm_svm_array(k,l);
    end
    mcc_den_2_part2=0;
    for f=1:1:length(cm_svm_array)
        if f ~=k
          for g=1:1:length(cm_svm_array)
            mcc_den_2_part2= mcc_den_2_part2+cm_svm_array(f,g);
          end
        end
    end
    mcc_denominator_2=(mcc_denominator_2+(mcc_den_2_part1*mcc_den_2_part2));
end

mcc = (mcc_numerator)/((mcc_denominator_1^0.5)*(mcc_denominator_2^0.5));