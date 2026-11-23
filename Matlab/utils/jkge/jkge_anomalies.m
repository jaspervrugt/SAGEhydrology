function [ay,aq] = jkge_anomalies(y,q,my,mq,good,method)
% Center long-term anomalies relative to an anchor to retain small
% variations.
ay=y(good)-my(good); aq=q(good)-mq(good);
if method==3
 fy=y(isfinite(y)); fq=q(isfinite(q));
 if ~isempty(fy),ay=(y(good)-fy(1))-mean(fy-fy(1));end
 if ~isempty(fq),aq=(q(good)-fq(1))-mean(fq-fq(1));end
end
end
