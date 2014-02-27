classdef Body331 < handle
    %KEVINTEST Summary of this class goes here
    %   Detailed explanation goes here
    properties
        Name
        Mass
        Radius
        % OrbitRadius and OrbitSpeed are given by means of both in real-life
        % perhaps will change to max radius and min speed in the future?
        OrbitRadius
        OrbitSpeed
        pos % [x;y]
        vel % [vx;vy]
        Color
        xHist
        yHist
        joined
    end
    methods
        function obj = Body331(pName,pOrbitRadius,pOrbitSpeed,pMass,pRadius,pColor)
            obj.Name = pName;
            obj.OrbitRadius = pOrbitRadius;
            obj.OrbitSpeed = pOrbitSpeed;
            randAngle = rand(1)*2*pi;
            obj.pos = [ obj.OrbitRadius * cos(randAngle); obj.OrbitRadius * sin(randAngle)];
            obj.vel = [ -1*obj.OrbitSpeed * cos(pi/2 - randAngle); obj.OrbitSpeed * sin(pi/2 - randAngle)];
            obj.Mass = pMass;
            obj.Radius = pRadius;
            obj.Color = pColor;
            obj.xHist = [];
            obj.yHist = [];
            obj.joined = false;
        end
        %returns a body for RK method.
        function newBody = rkCopy(obj,dx)
%             newBody = Body( obj.name, ...
%                             obj.OrbitRadius, ...
%                             obj.OrbitSpeed, ...
%                             obj.Mass, ...
%                             obj.Radius, ...
%                             obj.Color ...
%                             );
            newBody = Body331('Temporary', 0,0,obj.Mass,obj.Radius,obj.Color);
            newBody.vel = obj.vel;
            newBody.pos = obj.pos + dx;
        end

    end
end