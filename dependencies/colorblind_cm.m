function cm=colorblind_cm()
%Cal for colormap() function similar to 'bone' or any other colormap
cm=flipud(colorblind(size(bone,1)));