classdef SO3_SE3
    % Functions for handling SO3 and SE3 Lie Groups
    
    methods(Static)
        
        function out = tilde(x) %Eq.(9)
            % does both tilde and reverse tilde
            if min(size(x)) == 1
                out = [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];
            else
                out = 0.5*[x(3,2)-x(2,3); x(1,3)-x(3,1); x(2,1)-x(1,2)];
            end
        end
        
        function [R] = quat2R(q)
            % unit quaternion to rotation matrix
            e0 = q(1);
            e = q(2:4);
            R = eye(3) + 2*e0*SO3_SE3.tilde(e) + 2*SO3_SE3.tilde(e)*SO3_SE3.tilde(e);
            R(abs(R)<=eps)=0;
        end
        
        function [q] = R2quat(R)
            % rotation matrix to unit quaternion
            rotation_par = 1/sqrt(trace(R)+1)*(R-R');
            p = SO3_SE3.tilde(rotation_par);
            e0 = sqrt(1-p'*p/4);
            e = p/2;
            q = [e0;e];
        end
        
        function [q] = quatmult(qa, qb)
            e0 = qa(1)*qb(1) - dot(qa(2:4),qb(2:4));
            e = qa(1)*qb(2:4) + qb(1)*qa(2:4) + SO3_SE3.tilde(qa(2:4))*qb(2:4);
            q = [e0;e];
        end
        
        function R = ExpSO3(Omega) %Eq.(A.5)
            nOmega = norm(Omega);
            if nOmega == 0
                R = eye(3);
            else
                alpha = sin(nOmega)/nOmega;
                beta = 2*(1-cos(nOmega))/nOmega^2;
                R = eye(3) + alpha*SO3_SE3.tilde(Omega) + 0.5*beta*SO3_SE3.tilde(Omega)^2;
            end
        end
        
        function omega = LogSO3(R) %Eq.(A.8)
            k = acos(0.5*(trace(R)-1));
            if k == 0
                omega = zeros(3,1);
            elseif k < 1e-6
                omega = SO3_SE3.tilde(0.5*(R-R'));
            else
                k = k/sin(k);
                omega = SO3_SE3.tilde((0.5*k)*(R-R'));
            end
        end
        
        function Tm1 = TSO3m1(Omega) %Eq.(A.7)
            nOmega = norm(Omega);
            if nOmega == 0
                Tm1 = eye(3);
            else
                alpha = sin(nOmega)/nOmega;
                beta = 2*(1-cos(nOmega))/nOmega^2;
                Tm1 = eye(3) + SO3_SE3.tilde(0.5*Omega) + (1-alpha/beta)*1/nOmega^2*SO3_SE3.tilde(Omega)^2;
            end
        end
        
        function H_inv = InvSE3(H)
            R = H(1:3,1:3); x = H(1:3,4);
            H_inv = [R' -R'*x; 0 0 0 1];
        end
        
        function h = LogSE3(H) %Eq.(A.15)
            h = zeros(6,1);
            R = H(1:3,1:3); x = H(1:3,4);
            h(4:6) = SO3_SE3.LogSO3(R);
            h(1:3) = SO3_SE3.TSO3m1(h(4:6))'*x;
            h(abs(h)<1e-15) = 0;
        end
        
        function p = get_parameters_from_frame(H)
            % Euler Rodrigues Motion Parametrisation
            p = zeros(6,1);
            R = H(1:3,1:3); x = H(1:3,4);
            q = SO3_SE3.R2quat(R);
            p(4:6) = 2. * q(2:4);
            Tm1T = q(1) * eye(3) - SO3_SE3.tilde(q(2:4));
            p(1:3) = Tm1T * x;
        end
    end
end

