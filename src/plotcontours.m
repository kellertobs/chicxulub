function plotcontours(x,z,unit,wat,air)

if ~any(isnan(unit(:)))
    for i = 1:size(unit,3)
        contour(x,z,unit(:,:,i),1,'w','LineWidth',0.5);
    end
end

if any(wat(:))
    contour(x,z,wat,1,'w--','LineWidth',0.5);
end
if any(air(:))
    contour(x,z,air,1,'w:' ,'LineWidth',0.5);
end

end