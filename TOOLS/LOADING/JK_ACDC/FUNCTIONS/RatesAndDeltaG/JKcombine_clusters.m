function cn=JKcombine_clusters(c1,c2)
s1=regexp(c1,'[\D]+','match');
n1=str2double(regexp(c1,'[\d.]+','match'));
s2=regexp(c2,'[\D]+','match');
n2=str2double(regexp(c2,'[\d.]+','match'));

cn="";
for i=unique([s1,s2])
    c=0;
    for mol=1:size(s1,2)
        if i==s1(mol)
            c=c+n1(mol);
        end
    end
    for mol=1:size(s2,2)
        if i==s2(mol)
            c=c+n2(mol);
        end
    end
    cn=cn+num2str(c)+i;
end
end 