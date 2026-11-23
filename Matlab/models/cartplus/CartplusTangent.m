% ========================================================================
% SAGEhydrology | Model source
% File: Matlab/models/cartplus/CartplusTangent.m
%
% CaRT / CaRT+ numerical model and supporting routines.
%
% SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
% Upstream model: Xia (2026). Preserve the adjacent license notices.
% ========================================================================

classdef CartplusTangent
%CARTPLUSTANGENT Pair model values with one forward-mode derivative.
% A scalar object carries arrays of values and matching tangents. Only
% operations used by the standalone CaRT+ model are overloaded.

    properties (SetAccess = private)
        value
        slope
    end

    methods
        function x = CartplusTangent(value,slope)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CARTPLUSTANGENT Construct a primal array with forward derivative
        % values.
        %
        % SYNOPSIS:
        % CartplusTangent(...)
        %
        % INPUT ARGUMENTS:
        % value : Primal scalar or numerical operand.
        % slope : Forward derivative or vector of parameter-direction
        % derivatives.
        %
        % OUTPUT ARGUMENTS:
        % x : Computed numerical result; tangent lanes are retained.
        %
        % NOTES:
        % CARTPLUSTANGENT Seed a direction or promote a numeric constant.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nargin < 2
                slope = zeros(size(value));
            end
            x.value = value;
            x.slope = slope;
        end

        function y = real(x)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % REAL Extract the primal values from a tangent-array container.
        %
        % SYNOPSIS:
        % real(...)
        %
        % INPUT ARGUMENTS:
        % x : Numerical operand; scalar, vector, or forward-mode scalar.
        %
        % OUTPUT ARGUMENTS:
        % y : Computed numerical result; tangent lanes are retained.
        %
        % NOTES:
        % REAL Select branches using the primal value.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            y = x.value;
        end

        function y = size(x,varargin)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SIZE Return the dimensions of the primal numerical array.
        %
        % SYNOPSIS:
        % size(...)
        %
        % INPUT ARGUMENTS:
        % x : Numerical operand; scalar, vector, or forward-mode scalar.
        % varargin : Additional MATLAB operands or dimension arguments.
        %
        % OUTPUT ARGUMENTS:
        % y : Computed numerical result; tangent lanes are retained.
        %
        % NOTES:
        % SIZE Report the shape of the underlying array.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            y = size(x.value,varargin{:});
        end

        function y = end(x,k,n)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % END Return an array dimension for MATLAB end-indexing.
        %
        % SYNOPSIS:
        % end(...)
        %
        % INPUT ARGUMENTS:
        % x : Numerical operand; scalar, vector, or forward-mode scalar.
        % k : State index or release coefficient in this helper.
        % n : Number of intervals or state entries in this helper.
        %
        % OUTPUT ARGUMENTS:
        % y : Computed numerical result; tangent lanes are retained.
        %
        % NOTES:
        % END Resolve the last index for the underlying array.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dims = size(x.value);
            dims(end+1:n) = 1;
            if k == n
                y = prod(dims(k:end));
            else
                y = dims(k);
            end
        end

        function y = subsref(x,s)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SUBSREF Select matching primal entries and derivative lanes.
        %
        % SYNOPSIS:
        % subsref(...)
        %
        % INPUT ARGUMENTS:
        % x : Numerical operand; scalar, vector, or forward-mode scalar.
        % s : S used by this numerical operator.
        %
        % OUTPUT ARGUMENTS:
        % y : Computed numerical result; tangent lanes are retained.
        %
        % NOTES:
        % SUBSREF Take the same slice of the primal and tangent arrays.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(s(1).type,'()')
                index = s(1).subs;
                y = CartplusTangent(x.value(index{:}),x.slope(index{:}));
                if numel(s) > 1
                    y = subsref(y,s(2:end));
                end
            else
                y = builtin('subsref',x,s);
            end
        end

        function x = subsasgn(x,s,y)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SUBSASGN Assign matching primal entries and derivative lanes.
        %
        % SYNOPSIS:
        % subsasgn(...)
        %
        % INPUT ARGUMENTS:
        % x : Numerical operand; scalar, vector, or forward-mode scalar.
        % s : S used by this numerical operator.
        % y : Numerical state/operand; augmented with sensitivities when
        % enabled.
        %
        % OUTPUT ARGUMENTS:
        % x : Computed numerical result; tangent lanes are retained.
        %
        % NOTES:
        % SUBSASGN Assign the same slice of the primal and tangent arrays.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            y = CartplusTangent.lift(y);
            index = s.subs;
            x.value(index{:}) = y.value;
            x.slope(index{:}) = y.slope;
        end

        function y = plus(a,b)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PLUS Add primal values and derivative lanes.
        %
        % SYNOPSIS:
        % plus(...)
        %
        % INPUT ARGUMENTS:
        % a : First numerical operand or elevation-band fractions.
        % b : Second numerical operand or rate coefficient.
        %
        % OUTPUT ARGUMENTS:
        % y : Computed numerical result; tangent lanes are retained.
        %
        % NOTES:
        % PLUS Add values and tangents.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [a,b] = CartplusTangent.pair(a,b);
            y = CartplusTangent(a.value+b.value,a.slope+b.slope);
        end

        function y = uminus(a)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % UMINUS Negate primal values and derivative lanes.
        %
        % SYNOPSIS:
        % uminus(...)
        %
        % INPUT ARGUMENTS:
        % a : First numerical operand or elevation-band fractions.
        %
        % OUTPUT ARGUMENTS:
        % y : Computed numerical result; tangent lanes are retained.
        %
        % NOTES:
        % UMINUS Negate values and tangents.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            y = CartplusTangent(-a.value,-a.slope);
        end

        function y = minus(a,b)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MINUS Subtract primal values and derivative lanes.
        %
        % SYNOPSIS:
        % minus(...)
        %
        % INPUT ARGUMENTS:
        % a : First numerical operand or elevation-band fractions.
        % b : Second numerical operand or rate coefficient.
        %
        % OUTPUT ARGUMENTS:
        % y : Computed numerical result; tangent lanes are retained.
        %
        % NOTES:
        % MINUS Subtract values and tangents.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            y = plus(a,uminus(CartplusTangent.lift(b)));
        end

        function y = times(a,b)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TIMES Multiply elementwise using the forward product rule.
        %
        % SYNOPSIS:
        % times(...)
        %
        % INPUT ARGUMENTS:
        % a : First numerical operand or elevation-band fractions.
        % b : Second numerical operand or rate coefficient.
        %
        % OUTPUT ARGUMENTS:
        % y : Computed numerical result; tangent lanes are retained.
        %
        % NOTES:
        % TIMES Apply the elementwise product rule.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [a,b] = CartplusTangent.pair(a,b);
            slope = a.slope.*b.value+a.value.*b.slope;
            y = CartplusTangent(a.value.*b.value,slope);
        end

        function y = mtimes(a,b)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MTIMES Multiply operands using the forward product rule.
        %
        % SYNOPSIS:
        % mtimes(...)
        %
        % INPUT ARGUMENTS:
        % a : First numerical operand or elevation-band fractions.
        % b : Second numerical operand or rate coefficient.
        %
        % OUTPUT ARGUMENTS:
        % y : Computed numerical result; tangent lanes are retained.
        %
        % NOTES:
        % MTIMES Apply the matrix product rule.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [a,b] = CartplusTangent.pair(a,b);
            slope = a.slope*b.value+a.value*b.slope;
            y = CartplusTangent(a.value*b.value,slope);
        end

        function y = rdivide(a,b)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % RDIVIDE Divide elementwise using the forward quotient rule.
        %
        % SYNOPSIS:
        % rdivide(...)
        %
        % INPUT ARGUMENTS:
        % a : First numerical operand or elevation-band fractions.
        % b : Second numerical operand or rate coefficient.
        %
        % OUTPUT ARGUMENTS:
        % y : Computed numerical result; tangent lanes are retained.
        %
        % NOTES:
        % RDIVIDE Apply the elementwise quotient rule.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [a,b] = CartplusTangent.pair(a,b);
            v = a.value./b.value;
            y = CartplusTangent(v,(a.slope-v.*b.slope)./b.value);
        end

        function y = mrdivide(a,b)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MRDIVIDE Apply division with forward derivative propagation.
        %
        % SYNOPSIS:
        % mrdivide(...)
        %
        % INPUT ARGUMENTS:
        % a : First numerical operand or elevation-band fractions.
        % b : Second numerical operand or rate coefficient.
        %
        % OUTPUT ARGUMENTS:
        % y : Computed numerical result; tangent lanes are retained.
        %
        % NOTES:
        % MRDIVIDE Divide by a scalar as in the model equations.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            y = rdivide(a,b);
        end

        function y = power(a,b)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % POWER Raise values to elementwise powers and differentiate.
        %
        % SYNOPSIS:
        % power(...)
        %
        % INPUT ARGUMENTS:
        % a : First numerical operand or elevation-band fractions.
        % b : Second numerical operand or rate coefficient.
        %
        % OUTPUT ARGUMENTS:
        % y : Computed numerical result; tangent lanes are retained.
        %
        % NOTES:
        % POWER Differentiate a nonnegative base and variable exponent.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [a,b] = CartplusTangent.pair(a,b);
            v = a.value.^b.value;
            base = a.value;
            base(base == 0) = 1;
            slope = v.*(b.value.*a.slope./base+b.slope.*log(base));
            y = CartplusTangent(v,slope);
        end

        function y = mpower(a,b)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MPOWER Apply the scalar power operation with forward derivatives.
        %
        % SYNOPSIS:
        % mpower(...)
        %
        % INPUT ARGUMENTS:
        % a : First numerical operand or elevation-band fractions.
        % b : Second numerical operand or rate coefficient.
        %
        % OUTPUT ARGUMENTS:
        % y : Computed numerical result; tangent lanes are retained.
        %
        % NOTES:
        % MPOWER Raise a scalar store or coefficient to a power.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            y = power(a,b);
        end

        function y = ctranspose(a)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CTRANSPOSE Transpose primal values and derivative lanes.
        %
        % SYNOPSIS:
        % ctranspose(...)
        %
        % INPUT ARGUMENTS:
        % a : First numerical operand or elevation-band fractions.
        %
        % OUTPUT ARGUMENTS:
        % y : Computed numerical result; tangent lanes are retained.
        %
        % NOTES:
        % CTRANSPOSE Transpose primal and tangent arrays.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            y = CartplusTangent(a.value',a.slope');
        end

        function y = horzcat(varargin)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % HORZCAT Concatenate primal arrays and derivatives horizontally.
        %
        % SYNOPSIS:
        % horzcat(...)
        %
        % INPUT ARGUMENTS:
        % varargin : Additional MATLAB operands or dimension arguments.
        %
        % OUTPUT ARGUMENTS:
        % y : Computed numerical result; tangent lanes are retained.
        %
        % NOTES:
        % HORZCAT Assemble model outputs and spin-up stores.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            args = cellfun(@CartplusTangent.lift,varargin, ...
                'UniformOutput',false);
            values = cellfun(@(x) x.value,args,'UniformOutput',false);
            slopes = cellfun(@(x) x.slope,args,'UniformOutput',false);
            y = CartplusTangent([values{:}],[slopes{:}]);
        end

        function y = vertcat(varargin)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % VERTCAT Concatenate primal arrays and derivatives vertically.
        %
        % SYNOPSIS:
        % vertcat(...)
        %
        % INPUT ARGUMENTS:
        % varargin : Additional MATLAB operands or dimension arguments.
        %
        % OUTPUT ARGUMENTS:
        % y : Computed numerical result; tangent lanes are retained.
        %
        % NOTES:
        % VERTCAT Assemble memory and routing vectors.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            args = cellfun(@CartplusTangent.lift,varargin, ...
                'UniformOutput',false);
            values = cellfun(@(x) x.value,args,'UniformOutput',false);
            slopes = cellfun(@(x) x.slope,args,'UniformOutput',false);
            y = CartplusTangent(vertcat(values{:}),vertcat(slopes{:}));
        end

        function y = exp(a)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % EXP Evaluate exp and propagate forward derivatives.
        %
        % SYNOPSIS:
        % exp(...)
        %
        % INPUT ARGUMENTS:
        % a : First numerical operand or elevation-band fractions.
        %
        % OUTPUT ARGUMENTS:
        % y : Computed numerical result; tangent lanes are retained.
        %
        % NOTES:
        % EXP Differentiate the exponential.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            v = exp(a.value);
            y = CartplusTangent(v,v.*a.slope);
        end

        function y = expm1(a)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % EXPM1 Evaluate expm1 and propagate forward derivatives.
        %
        % SYNOPSIS:
        % expm1(...)
        %
        % INPUT ARGUMENTS:
        % a : First numerical operand or elevation-band fractions.
        %
        % OUTPUT ARGUMENTS:
        % y : Computed numerical result; tangent lanes are retained.
        %
        % NOTES:
        % EXPM1 Differentiate the accurate small-argument exponential.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            y = CartplusTangent(expm1(a.value),exp(a.value).*a.slope);
        end

        function y = log(a)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % LOG Evaluate log and propagate forward derivatives.
        %
        % SYNOPSIS:
        % log(...)
        %
        % INPUT ARGUMENTS:
        % a : First numerical operand or elevation-band fractions.
        %
        % OUTPUT ARGUMENTS:
        % y : Computed numerical result; tangent lanes are retained.
        %
        % NOTES:
        % LOG Differentiate the natural logarithm.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            y = CartplusTangent(log(a.value),a.slope./a.value);
        end

        function y = log1p(a)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % LOG1P Evaluate log1p and propagate forward derivatives.
        %
        % SYNOPSIS:
        % log1p(...)
        %
        % INPUT ARGUMENTS:
        % a : First numerical operand or elevation-band fractions.
        %
        % OUTPUT ARGUMENTS:
        % y : Computed numerical result; tangent lanes are retained.
        %
        % NOTES:
        % LOG1P Differentiate the accurate small-argument logarithm.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            y = CartplusTangent(log1p(a.value),a.slope./(1+a.value));
        end

        function y = tanh(a)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TANH Evaluate tanh and propagate forward derivatives.
        %
        % SYNOPSIS:
        % tanh(...)
        %
        % INPUT ARGUMENTS:
        % a : First numerical operand or elevation-band fractions.
        %
        % OUTPUT ARGUMENTS:
        % y : Computed numerical result; tangent lanes are retained.
        %
        % NOTES:
        % TANH Differentiate the hyperbolic tangent.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            v = tanh(a.value);
            y = CartplusTangent(v,(1-v.^2).*a.slope);
        end

        function y = mean(a)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MEAN Average primal values and their derivative lanes.
        %
        % SYNOPSIS:
        % mean(...)
        %
        % INPUT ARGUMENTS:
        % a : First numerical operand or elevation-band fractions.
        %
        % OUTPUT ARGUMENTS:
        % y : Computed numerical result; tangent lanes are retained.
        %
        % NOTES:
        % MEAN Average primal values and slopes along the same dimension.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            y = CartplusTangent(mean(a.value),mean(a.slope));
        end

        function y = diff(a)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DIFF Difference consecutive primal values and derivative lanes.
        %
        % SYNOPSIS:
        % diff(...)
        %
        % INPUT ARGUMENTS:
        % a : First numerical operand or elevation-band fractions.
        %
        % OUTPUT ARGUMENTS:
        % y : Computed numerical result; tangent lanes are retained.
        %
        % NOTES:
        % DIFF Difference adjacent primal values and slopes.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            y = CartplusTangent(diff(a.value),diff(a.slope));
        end

        function y = conv(a,b)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CONV Convolve values while applying the forward product rule.
        %
        % SYNOPSIS:
        % conv(...)
        %
        % INPUT ARGUMENTS:
        % a : First numerical operand or elevation-band fractions.
        % b : Second numerical operand or rate coefficient.
        %
        % OUTPUT ARGUMENTS:
        % y : Computed numerical result; tangent lanes are retained.
        %
        % NOTES:
        % CONV Apply the product rule to unit-hydrograph routing.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [a,b] = CartplusTangent.pair(a,b);
            slope = conv(a.slope,b.value)+conv(a.value,b.slope);
            y = CartplusTangent(conv(a.value,b.value),slope);
        end
    end

    methods (Static)
        function x = lift(x)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % LIFT Promote an ordinary value to a zero-derivative tangent.
        %
        % SYNOPSIS:
        % lift(...)
        %
        % INPUT ARGUMENTS:
        % x : Numerical operand; scalar, vector, or forward-mode scalar.
        %
        % OUTPUT ARGUMENTS:
        % x : Computed numerical result; tangent lanes are retained.
        %
        % NOTES:
        % LIFT Promote numeric constants to zero-tangent values.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~isa(x,'CartplusTangent')
                x = CartplusTangent(x);
            end
        end

        function [a,b] = pair(a,b)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PAIR Promote two operands to compatible tangent-array shapes.
        %
        % SYNOPSIS:
        % pair(...)
        %
        % INPUT ARGUMENTS:
        % a : First numerical operand or elevation-band fractions.
        % b : Second numerical operand or rate coefficient.
        %
        % OUTPUT ARGUMENTS:
        % a : First promoted numerical operand.
        % b : Second promoted numerical operand.
        %
        % NOTES:
        % PAIR Promote both operands before applying a derivative rule.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            a = CartplusTangent.lift(a);
            b = CartplusTangent.lift(b);
        end
    end
end
