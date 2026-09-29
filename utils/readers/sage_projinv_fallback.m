function [latitude,longitude] = sage_projinv_fallback(epsg,x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SAGE_PROJINV_FALLBACK Convert supported projected basin coordinates.
%
% SYNOPSIS:
%   [latitude,longitude] = sage_projinv_fallback(epsg,x,y)
%
% INPUT ARGUMENTS:
%   epsg            projected CRS: 3035, 3857, 25832, or 3057
%   x               projected easting in meters
%   y               projected northing in meters
%
% OUTPUT ARGUMENTS:
%   latitude        geographic latitude in degrees
%   longitude       geographic longitude in degrees
%
% NOTES:
%   Uses the inverse equations for the four regional SAGE schemas.
%   This is not a general-purpose replacement for Mapping Toolbox.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~isequal(size(x),size(y))
        error('SAGE:projectionSizeMismatch', ...
            ['Projected x and y coordinates ' ...
            'must have the same size.']);
    end
    x = double(x);
    y = double(y);
    a = 6378137;
    f = 1 / 298.257222101;  % GRS80: EPSG 3035, 25832, 3057
    e2 = f * (2 - f);
    e = sqrt(e2);

    switch double(epsg)
        case 3857             % WGS 84 / Pseudo-Mercator
            longitude = (x / a) * 180 / pi;
            latitude = atan(sinh(y / a)) * 180 / pi;

        case 3035             % ETRS89-extended / LAEA Europe
            phi0 = 52 * pi / 180;
            lambda0 = 10 * pi / 180;
            qp = local_authalic_q(pi / 2,e2,e);
            beta0 = asin(local_authalic_q(phi0,e2,e) / qp);
            rq = a * sqrt(qp / 2);
            m0 = cos(phi0) / sqrt(1 - e2 * sin(phi0)^2);
            D = a * m0 / (rq * cos(beta0));
            xx = (x - 4321000) / D;
            yy = (y - 3210000) * D;
            rho = hypot(xx,yy);
            c = 2 * asin(min(rho / (2 * rq),1));
            safeRho = rho;
            safeRho(rho == 0) = 1;
            beta = asin(max(-1,min(1, ...
                cos(c) * sin(beta0) ...
                + yy .* sin(c) * cos(beta0) ./ safeRho)));
            beta(rho == 0) = beta0;
            lambda = lambda0 + atan2(xx .* sin(c), ...
                rho * cos(beta0) .* cos(c) ...
                - yy * sin(beta0) .* sin(c));
            lambda(rho == 0) = lambda0;
            phi = local_inverse_authalic_q(qp * sin(beta),e2,e);
            latitude = phi * 180 / pi;
            longitude = lambda * 180 / pi;

        case 25832            % ETRS89 / UTM zone 32N
            k0 = 0.9996;
            ep2 = e2 / (1 - e2);
            M = y / k0;
            mu = M / (a * (1 - e2/4 - 3*e2^2/64 ...
                - 5*e2^3/256));
            e1 = (1 - sqrt(1 - e2)) / (1 + sqrt(1 - e2));
            phi1 = mu ...
                + (3*e1/2 - 27*e1^3/32) .* sin(2*mu) ...
                + (21*e1^2/16 - 55*e1^4/32) .* sin(4*mu) ...
                + (151*e1^3/96) .* sin(6*mu) ...
                + 1097*e1^4/512 .* sin(8*mu);
            sinPhi1 = sin(phi1);
            cosPhi1 = cos(phi1);
            tanPhi1 = tan(phi1);
            N1 = a ./ sqrt(1 - e2 * sinPhi1.^2);
            R1 = a * (1 - e2) ./ ...
                (1 - e2 * sinPhi1.^2).^(3/2);
            T1 = tanPhi1.^2;
            C1 = ep2 * cosPhi1.^2;
            D = (x - 500000) ./ (N1 * k0);
            phi = phi1 - (N1 .* tanPhi1 ./ R1) .* ...
                (D.^2/2 ...
                - (5 + 3*T1 + 10*C1 - 4*C1.^2 - 9*ep2) ...
                .* D.^4/24 ...
                + (61 + 90*T1 + 298*C1 + 45*T1.^2 ...
                - 252*ep2 - 3*C1.^2) .* D.^6/720);
            lambda = 9*pi/180 + ...
                (D - (1 + 2*T1 + C1) .* D.^3/6 ...
                + (5 - 2*C1 + 28*T1 - 3*C1.^2 ...
                + 8*ep2 + 24*T1.^2) .* D.^5/120) ...
                ./ cosPhi1;
            latitude = phi * 180 / pi;
            longitude = lambda * 180 / pi;

        case 3057             % ISN93 / Lambert 1993
            phi0 = 65*pi/180;
            phi1 = 64.25*pi/180;
            phi2 = 65.75*pi/180;
            lambda0 = -19*pi/180;
            m1 = local_lcc_m(phi1,e2);
            m2 = local_lcc_m(phi2,e2);
            t1 = local_lcc_t(phi1,e);
            t2 = local_lcc_t(phi2,e);
            n = (log(m1) - log(m2)) / (log(t1) - log(t2));
            F = m1 / (n * t1^n);
            rho0 = a * F * local_lcc_t(phi0,e)^n;
            xx = x - 500000;
            yy = rho0 - (y - 500000);
            rho = hypot(xx,yy);
            t = (rho / (a * F)).^(1/n);
            phi = pi/2 - 2*atan(t);
            for k = 1:8
                phi = pi/2 - 2*atan(t .* ...
                    ((1 - e*sin(phi)) ./ ...
                    (1 + e*sin(phi))).^(e/2));
            end
            lambda = lambda0 + atan2(xx,yy) / n;
            latitude = phi * 180 / pi;
            longitude = lambda * 180 / pi;

        otherwise
            error('SAGE:unsupportedProjection', ...
                ['No Mapping Toolbox fallback ' ...
                'is defined for EPSG:%g.'], ...
                double(epsg));
    end
end

function q = local_authalic_q(phi,e2,e)
    s = sin(phi);
    q = (1 - e2) * (s ./ (1 - e2*s.^2) ...
        - log((1 - e*s) ./ (1 + e*s)) / (2*e));
end

function phi = local_inverse_authalic_q(q,e2,e)
    qp = local_authalic_q(pi/2,e2,e);
    phi = asin(max(-1,min(1,q / qp)));
    for k = 1:8
        dq = 2 * (1 - e2) * cos(phi) ./ ...
            (1 - e2*sin(phi).^2).^2;
        phi = phi + (q - local_authalic_q(phi,e2,e)) ./ dq;
    end
end

function m = local_lcc_m(phi,e2)
    m = cos(phi) ./ sqrt(1 - e2*sin(phi).^2);
end

function t = local_lcc_t(phi,e)
    t = tan(pi/4 - phi/2) ./ ...
        ((1 - e*sin(phi)) ./ (1 + e*sin(phi))).^(e/2);
end
