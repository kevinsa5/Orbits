classdef Body < handle
    properties
        Name
        Mass
        Radius
        OrbitRadius
        OrbitSpeed
        pos % [x;y;z]
        vel % [vx;vy;vz]
        Color
        xHist
        yHist
        joined %crashed into someone else
        originalName
    end
    methods
        function obj = Body(pName,pPos,pVel,pMass,pRadius,pColor)
            obj.Name = pName;
            obj.pos = pPos;
            obj.vel = pVel;
            obj.Mass = pMass;
            obj.Radius = pRadius;
            obj.Color = pColor;
            obj.xHist = [];
            obj.yHist = [];
            obj.joined = false;
            obj.originalName = pName;
        end
%         function moon = makeMoon(obj,pName, pOrbitRadius, pOrbitSpeed, pMass, pRadius, pColor)
%             moon = Body(pName,pOrbitRadius, pOrbitSpeed, pMass, pRadius, pColor);
%             randAngle = rand(1)*2*pi;
%             moon.pos = obj.pos + [moon.OrbitRadius * cos(randAngle); moon.OrbitRadius * sin(randAngle)];
%             moon.vel = obj.vel + [-1*moon.OrbitSpeed * cos(pi/2 - randAngle); moon.OrbitSpeed * sin(pi/2 - randAngle)];
%             obj.numberOfChildren = obj.numberOfChildren + 1;
%             moon.childNumber = obj.numberOfChildren;
%             moon.parent = obj;
%         end
        %returns a body for RK method.
        function newBody = rkCopy(obj,dx)
            newBody = Body('Temporary', obj.pos,obj.vel,obj.Mass,obj.Radius,obj.Color);
            newBody.pos = obj.pos + dx;
        end

    end
end