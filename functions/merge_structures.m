function finalstruct=merge_structures(structure1,structure2)

% function to merge two structures together
finalstruct=structure1;
f=fieldnames(structure2);

for I=1:length(f)
    finalstruct.(f{I})=structure2.(f{I});
end