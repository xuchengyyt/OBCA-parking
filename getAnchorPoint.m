function[list] = getAnchorPoint(Deriction)
list = [];
for i = 1:length(Deriction)-1
    if(Deriction(i+1) * Deriction(i)<0)
        list = [list;i];
    end
end
end