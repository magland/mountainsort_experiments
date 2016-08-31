function compile_results

[alg_names,ds_names]=get_alg_ds_names
unit_numbers.h1=[1];
unit_numbers.m1=[];
unit_numbers.m2=[];
unit_numbers.m3=[];
unit_numbers.m4=[];
unit_numbers.m5=[];

alg_names={'ms80'};
ds_names={'m1'};

for a=1:length(alg_names)
    for b=1:length(ds_names)
        algname=alg_names{a};
        dsname=ds_names{b};
        output_path=['output/',algname,'-',dsname];
        CM=readmda([output_path,'/confusion_matrix.csv']);
        LM=readmda([output_path,'/optimal_label_map.csv']);
        stats=compute_stats(CM,LM);
        tmp=unit_numbers.(dsname);
        for ii=1:length(tmp)
            k=tmp(ii);
            disp(CM);
            disp(stats(1:3,:));
            LM
            fprintf('%s %s %d num=%g, fn=%g, fp=%g\n',algname,dsname,k,stats(1,k),stats(4,k),stats(4,k));
        end;
    end;
end;

end

function [alg_names,ds_names]=get_alg_ds_names

alg_names={};
ds_names={};

tmp=dir('output');
for j=1:length(tmp)
    folder=tmp(j);
    if (folder.isdir)
        if (folder.name(1)~='.')
            vals=strsplit(folder.name,'-');
            if (~contains_string(alg_names,vals{1}))
                alg_names{end+1}=vals{1};
            end;
            if (~contains_string(ds_names,vals{2}))
                ds_names{end+1}=vals{2};
            end;
        end;
    end;
end;

end

function ret=contains_string(names,str)
ret=0;
for j=1:length(names)
    if (strcmp(names{j},str))
        ret=1;
    end;
end;
end

function ret=compute_stats(confusion_matrix,label_map)

K1=size(confusion_matrix,1)-1;
K2=size(confusion_matrix,2)-1;

inv_label_map=zeros(1,K1);
for k=1:K2
    if (label_map(k)) inv_label_map(label_map(k))=k; end;
end;

%total, false negatives, false positives
ret=zeros(5,K1);
for k1=1:K1
    if (label_map(k1))
        a=sum(confusion_matrix(k1,:));
        b=confusion_matrix(k1,inv_label_map(k1));
        c=sum(confusion_matrix(:,inv_label_map(k1)));
        ret(1,k1)=a; %total
        ret(2,k1)=a-b; %false negatives
        ret(3,k1)=c-b; %false positives
    else
        a=sum(confusion_matrix(k1,:));
        ret(1,k1)=a; %total
        ret(2,k1)=a; %false negatives
        ret(3,k1)=0; %false positives
    end;
    if (ret(1,k1))
        ret(4,k1)=ret(2,k1)/ret(1,k1);
        ret(5,k1)=ret(3,k1)/ret(1,k1);
    end;
end;

end
