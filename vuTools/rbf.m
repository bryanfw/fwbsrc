function weight = rbf(x,xc,scale)
    r = sqrt(sum((x-repmat(xc,length(x),1)).^2,2))./scale;
    weight = sparse(zeros(length(x),1));
    idx = find((r<repmat(1,length(x),1))==1);
    weight(idx) = (1-r(idx)).^4.*(3.*r(idx).^3 + 12.*r(idx).^2 + 16.*r(idx) + 4)./4.;    
end