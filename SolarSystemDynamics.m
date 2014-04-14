function varargout = SolarSystemDynamics(varargin)
    % Edit the above text to modify the response to help SolarSystemDynamics

    % Last Modified by GUIDE v2.5 11-Apr-2014 19:13:59

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @SolarSystemDynamics_OpeningFcn, ...
                       'gui_OutputFcn',  @SolarSystemDynamics_OutputFcn, ...
                       'gui_LayoutFcn',  [] , ...
                       'gui_Callback',   []);
    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1});
    end

    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end
    % End initialization code - DO NOT EDIT
end

% --- Executes just before SolarSystemDynamics is made visible.
function SolarSystemDynamics_OpeningFcn(hObject, ~, handles, varargin)
    % Choose default command line output for SolarSystemDynamics
    handles.output = hObject;
    % Update handles structure
    guidata(hObject, handles);
    % UIWAIT makes SolarSystemDynamics wait for user response (see UIRESUME)
    % uiwait(handles.MainFigure);
    reset(handles);
%    clearvars;
    clc;
%    clear all;
    cla;
    global bodies;
    bodies = [];
    
    menuConfiguration_Callback(0,0,handles);
    menuLaunchDate_Callback(handles.menuLaunchDate, 0, handles)
    h = zoom;
    set(h,'ActionPostCallback',@drawLabels);
    
end

% --- Outputs from this function are returned to the command line.
function varargout = SolarSystemDynamics_OutputFcn(~, ~, handles) 
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % ~  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.output;
end

function drawBodies(handles)
    global bodies
    global spaceship;
    %global planetLabeling; planetLabeling = get(handles.chkLabeling, 'Value');
    %global moonLabeling; moonLabeling = get(handles.chkMoonLabeling, 'Value');

    lims = axis;
    contents = cellstr(get(handles.centeredBody,'String'));
    centeredBodyName = contents{get(handles.centeredBody,'Value')};
    
    translation = [0;0];
    for i = bodies
        if strcmp(centeredBodyName, i.originalName), translation = i.pos; end
    end
    
    cla;
    plot(spaceship.xHist,spaceship.yHist);
    hold on;
    
    for i=1:length(bodies)
        if bodies(i).joined, continue; end
        % using rectangle allows doing circle diameter in AU instead of px
        corner = transpose(bodies(i).pos - bodies(i).Radius);
        %disp(bodies(i).Name);
        rectangle('Position',[corner(1:2) 2*bodies(i).Radius 2*bodies(i).Radius],'Curvature',1);
        plot(bodies(i).xHist, bodies(i).yHist, bodies(i).Color);
    end
    % used for centering on the currently tracked body:
    deltaAxis = [translation(1) - mean(lims(1:2)), translation(2) - mean(lims(3:4))];
    axis(lims + [deltaAxis(1) deltaAxis(1) deltaAxis(2) deltaAxis(2)]);
    
    hold off;
end

function drawLabels(~,~)
    handles = guidata(gcf);
    planetLabeling = get(handles.chkLabeling, 'Value');
    moonLabeling = get(handles.chkMoonLabeling, 'Value');
    
    if strcmp(get(handles.btnGo, 'String'), 'Go') || ...
       strcmp(get(handles.btnGo, 'String'), 'Resume')
        drawBodies(handles);
    end
    global bodies;
    global planetList;
    global moonList;
    
    if planetLabeling || moonLabeling
        yax = ylim;
        xax = xlim;
        dx = 0.02*(xax(2)-xax(1));
        dy = 0.02*(yax(2)-yax(1));
        for i=1:length(bodies)
            if bodies(i).joined, continue; end
            if planetLabeling && sum(sum(strcmp(bodies(i).Name,planetList)))
                text(bodies(i).pos(1)+dx, bodies(i).pos(2)+dy, bodies(i).Name);
            end
            if moonLabeling && sum(sum(strcmp(bodies(i).Name,moonList)))
                text(bodies(i).pos(1)+dx, bodies(i).pos(2)+dy, bodies(i).Name);
            end         
        end
        t = str2double(get(handles.FrameCount, 'String'));
        days = 365.242*t*str2double(get(handles.txtTimeStep, 'String'));
        text(xax(1)+dx, yax(2)-dy, strcat(num2str(days),' days'));
    end
end

function btnGo_Callback(hObject, ~, handles) %#ok<DEFNU>
    global bodies;
    global spaceship;
    global venus;
    capturingMovie = false;
    if capturingMovie
        vidWriter = VideoWriter('movies/movie.avi');
        open(vidWriter);
    end

    if strcmp(get(hObject, 'String'), 'Pause')
        disp('pause was hit');
        set(hObject, 'String', 'Resume');
        set(handles.btnReset, 'Enable','on');
        return;
    end

    disp('started');
    start = str2double(get(hObject,'UserData'));

    set(hObject, 'String', 'Pause');
    set(handles.menuConfiguration, 'Enable', 'off');
    set(handles.txtTimeStep, 'Enable','off');
    set(handles.txtFrameSkip, 'Enable','off');
    set(handles.btnReset, 'Enable','off');
    set(handles.menuMethod, 'Enable', 'off');
    set(handles.menuLaunchDate, 'Enable', 'off');

    deltaT = str2double(get(handles.txtTimeStep, 'String'));
    frameSkip = str2double(get(handles.txtFrameSkip, 'String'));

    % 100000 is the arbitrary number of iterations I picked
    for t=start:10000000
        % window was closed?
        if ~ ishandle(handles.MainFigure)
            break;
        end

        if strcmp(get(hObject, 'String'), 'Resume')
            set(hObject, 'UserData', num2str(t));
            %drawBodies(true);
            drawBodies(handles);
            drawLabels(0,0);

            %axis(defaultAxes);
            return;
        end
        set(hObject, 'UserData', strcat('Running ', num2str(t)));
        set(handles.FrameCount, 'String', num2str(t));

        contents = cellstr(get(handles.menuMethod,'String'));
        method = contents{get(handles.menuMethod,'Value')};
        
        if strcmp(method,'Forward Euler')
            %eventually, do j=(i+1):length(bodies) and accelerate j based on i as well.
            %will have to change accel to force and back again
            %right now, we're duplicating all Gm/r^2 calculations
            for i=1:length(bodies)
                for j=1:length(bodies)
                    if i == j || bodies(i).joined || bodies(j).joined
                        continue
                    end
                    accel = forceOn(bodies(i),bodies(j));
                    bodies(i).vel = bodies(i).vel + accel*deltaT;
                end
            end
            % now, move all the bodies
            % can't combine with the acceleration loop because you can't calculate
            % forces on one, then move it, then use that new position for other
            % bodies' forces. move all, then accelerate all.
            for i=1:length(bodies)
                if bodies(i).joined, continue; end
                bodies(i).pos = bodies(i).pos + bodies(i).vel * deltaT;
                bodies(i).xHist = [bodies(i).xHist bodies(i).pos(1)];
                bodies(i).yHist = [bodies(i).yHist bodies(i).pos(2)];
            end
        elseif strcmp(method,'Runge Kutta 4')
            % see http://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#The_Runge.E2.80.93Kutta_method
            % This implementation is described in http://spiff.rit.edu/richmond/nbody/OrbitRungeKutta4.pdf

            %first, calculate k1 for all particles, then k2 for all particles,
            %etc
            order = 4;
            kr = zeros(3,order,length(bodies));
            kv = zeros(3,order,length(bodies));
            %set up kr1 and kv1 for all bodies:
            for i = 1:length(bodies)
                if bodies(i).joined, continue; end
                kr(:,1,i) = bodies(i).vel;
                accel = 0;
                for j = 1:length(bodies)
                    if bodies(j).joined, continue; end
                    if i ~= j 
                        accel = accel + forceOn(bodies(i),bodies(j));
                    end
                end
                kv(:,1,i) = accel;
            end
            % calculate kr2, kv2 for each
            for i = 1:length(bodies)
                if bodies(i).joined, continue; end
                kr(:,2,i) = bodies(i).vel + kv(:,1,i)*deltaT/2;
                accel = 0;
                for j = 1:length(bodies)
                    if bodies(j).joined, continue; end
                    if i ~= j
                        accel = accel + forceOn( bodies(i).rkCopy(deltaT/2 * kr(:,1,i)), ...
                                                 bodies(j).rkCopy(deltaT/2 * kr(:,1,j)));
                    end
                end
                kv(:,2,i) = accel;
            end
            % calculate kr3, kv3 for each 
            for i = 1:length(bodies)
                if bodies(i).joined, continue; end
                kr(:,3,i) = bodies(i).vel + kv(:,2,i)*deltaT/2;
                accel = 0;
                for j = 1:length(bodies)
                    if bodies(j).joined, continue; end
                    if i ~= j
                        accel = accel + forceOn( bodies(i).rkCopy(deltaT/2 * kr(:,2,i)), ...
                                                 bodies(j).rkCopy(deltaT/2 * kr(:,2,j)));
                    end
                end
                kv(:,3,i) = accel;
            end
            % calculate kr4, kv4 for each 
            for i = 1:length(bodies)
                if bodies(i).joined, continue; end
                kr(:,4,i) = bodies(i).vel + kv(:,3,i)*deltaT;
                accel = 0;
                for j = 1:length(bodies)
                    if bodies(j).joined, continue; end
                    if i ~= j
                        accel = accel + forceOn( bodies(i).rkCopy(deltaT * kr(:,3,i)), ...
                                                 bodies(j).rkCopy(deltaT * kr(:,3,j)));
                    end
                end
                kv(:,4,i) = accel;
            end
            for i=1:length(bodies)
                if bodies(i).joined, continue; end
                bodies(i).vel = bodies(i).vel + (deltaT/6) * (kv(:,1,i) + 2*kv(:,2,i) + 2*kv(:,3,i) + kv(:,4,i));
                bodies(i).pos = bodies(i).pos + (deltaT/6) * (kr(:,1,i) + 2*kr(:,2,i) + 2*kr(:,3,i) + kr(:,4,i));
            end
            for i=1:length(bodies)
                if bodies(i).joined, continue; end
                bodies(i).xHist = [bodies(i).xHist bodies(i).pos(1)];
                bodies(i).yHist = [bodies(i).yHist bodies(i).pos(2)];
            end
        elseif strcmp(method, 'Verlet')
            % Taken from http://xbeams.chem.yale.edu/~batista/vaa/node60.html
            dX = zeros(3,length(bodies));
            dV = zeros(3,length(bodies));
            for i = 1:length(bodies)
                if bodies(i).joined, continue; end
                accel = [0;0;0];
                iAccel = [0;0;0];
                for j = 1:length(bodies)
                    if bodies(j).joined, continue; end
                    if i ~= j
                        accel = accel + forceOn( bodies(i), bodies(j));
                        iAccel = iAccel + forceOn( bodies(i).rkCopy(deltaT*bodies(i).vel),...
                                                   bodies(j).rkCopy(deltaT*bodies(j).vel));
                    end
                end
                dX(:,i) = bodies(i).vel*deltaT + 0.5 * deltaT^2 * accel;
                dV(:,i) = 0.5 * (accel + iAccel) * deltaT;
            end

            for i=1:length(bodies)
                if bodies(i).joined, continue; end
                bodies(i).pos = bodies(i).pos + dX(:,i);
                bodies(i).vel = bodies(i).vel + dV(:,i);

                bodies(i).xHist = [bodies(i).xHist bodies(i).pos(1)];
                bodies(i).yHist = [bodies(i).yHist bodies(i).pos(2)];
            end
        end
        
        % check for collisions
        for i = 1:length(bodies)
            if bodies(i).joined, continue; end
           for j = 1:length(bodies)
               if bodies(j).joined, continue; end
               if i == j 
                   continue;
               end
               d = norm(bodies(i).pos - bodies(j).pos);
               if d < bodies(i).Radius + bodies(j).Radius
                   totalMomentum = bodies(i).Mass .* bodies(i).vel + bodies(j).Mass .* bodies(j).vel;
                   biggerIndex = i;
                   smallerIndex = j;
                   if bodies(j).Mass > bodies(i).Mass
                       biggerIndex = j;
                       smallerIndex = i;
                   end
                   bodies(biggerIndex).Name = strcat(bodies(biggerIndex).Name, '+', bodies(smallerIndex).Name);
                   bodies(biggerIndex).Mass = bodies(i).Mass + bodies(j).Mass;
                   bodies(biggerIndex).vel = totalMomentum / bodies(biggerIndex).Mass;
                   bodies(smallerIndex).joined = true;
               end
           end
        end
        
        if norm(spaceship.pos - venus.pos) < 0.01
            disp(norm(spaceship.pos - venus.pos));
        end

        %axis(defaultAxes);
        if mod(t,frameSkip) == 0 
            drawBodies(handles);
            drawLabels(0,0);
            drawnow;
            if capturingMovie && ishandle(handles.MainFigure)
                writeVideo(vidWriter,getframe(gca));
            end
        end
    end
    if ishandle(hObject)
        set(hObject, 'UserData', 'Finished');
        set(hObject, 'String', 'Go');
        set(hObject, 'UserData', '1');
    end
    
    if capturingMovie
        close(vidWriter);
    end
end

function accel = forceOn(body1, body2)
    % Big G in units of Au^3 / (earth mass * year^2)
    % extra sig figs, just for funsies
    G = 1.18556068632395e-04;
    r = norm(body1.pos - body2.pos);
    % next line has (x1-x2)/mag(x1-x2) to represent r-hat
    accel = (-G*body2.Mass * (body1.pos - body2.pos)) / r^3;
end
    
function txtTimeStep_Callback(hObject, ~, ~) %#ok<DEFNU>
    val = get(hObject, 'String');
    if isnan(str2double(val))
        set(hObject, 'String', 0.01);
        msgbox(strcat('Cannot convert ''', val, ''' to a number'),'Time Step not valid');
    elseif str2double(val) <= 0
        set(hObject, 'String', 0.01);
        msgbox('Time Step must be greater than zero', 'Time Step not valid');
    end
end

function txtFrameSkip_Callback(hObject, ~, ~) %#ok<DEFNU>
    val = get(hObject, 'String');
    if isnan(str2double(val))
        set(hObject, 'String', 1);
        msgbox(strcat('Cannot convert ''', val, ''' to a number'), 'Frame Skip not valid');
    elseif str2double(val) <= 0
        set(hObject, 'String', 1);
        msgbox('Frame Skip must be greater than zero', 'Frame Skip not valid');
    end
end

function MainFigure_CloseRequestFcn(hObject, ~, ~) %#ok<DEFNU>
    delete(hObject);
end

function reset(handles)
    global defaultAxes;

    defaultAxes = [-30 30 -30 30];
    axis(defaultAxes);
    set(handles.txtTimeStep, 'String', 0.000273);
    set(handles.txtFrameSkip, 'String', 1);
    set(handles.btnGo, 'UserData', '1');
    set(handles.FrameCount,'String',0);
    
    set(handles.menuConfiguration, 'Enable', 'on');
    set(handles.txtTimeStep, 'Enable','on');
    set(handles.txtFrameSkip, 'Enable','on');
    set(handles.menuMethod, 'Enable', 'on');
    set(handles.menuLaunchDate, 'Enable', 'on');
    set(handles.btnGo, 'String', 'Go');
end

function btnReset_Callback(~, ~, handles) %#ok<DEFNU>
    reset(handles);
end

function menuMethod_Callback(~, ~, ~) %#ok<DEFNU>
    % maybe set a good default timestep here?
end

function menuMethod_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

function chkLabeling_Callback(~, ~, handles) %#ok<DEFNU>
    if strcmp(get(handles.btnGo,'String'), 'Go') || ...
       strcmp(get(handles.btnGo,'String'), 'Resume')
        drawBodies(handles);
        drawLabels(0,0);
    end
end

function menuConfiguration_Callback(~, ~, handles) 
    global bodyData;
    global planetList;
    global moonList;
    contents = cellstr(get(handles.menuConfiguration,'String'));
    configuration = contents{get(handles.menuConfiguration,'Value')};
    
    %bodyPositions is a list of tables
    primaryData = readtable('body_data/info.txt');
    secondaryData = readtable('body_data/mooninfo_1.txt');
    planetList = primaryData.Name;
    moonList = primaryData.Name;
    % bodyInfo.txt expects data of the form of
    % BodyName, BodyMass, BodyRadius
    
    if strcmp(configuration, 'Full Solar System')
        bodyData = [primaryData; secondaryData];
    elseif strcmp(configuration, 'Sun, Planets only')
        disp('planets only!');
        bodyData = primaryData;
    elseif strcmp(configuration, 'Sun, Earth, Moon')
        bodyData = [primaryData; secondaryData];
        criteria = strcmp(bodyData.Name,'Sun') | strcmp(bodyData.Name,'Earth') | strcmp(bodyData.Name,'Moon');
        bodyData = bodyData(criteria,:);
    end
    for row = 1:size(bodyData,1)
        mass = bodyData.Mass(row); % e20 kg
        radius = bodyData.Radius(row); % km
        bodyData.Mass(row) = mass * 1.6744252e-5; %convert to earth masses
        bodyData.Radius(row) = radius * 6.685e-9;   %convert to au
    end
    % now, bodyData has only the data for the bodies we want
    nameList{1} = 'Ship';
    for i = 1:length(bodyData.Name)
        nameList{i+1} = bodyData.Name{i};
    end
    
    set(handles.centeredBody,'String',nameList);
    set(handles.centeredBody,'Value',1);
    
    sunFile = readtable('body_data/Sun.txt');
    
    set(handles.menuLaunchDate,'Enable','on');
    set(handles.menuLaunchDate,'String',sunFile.epoch);
    set(handles.menuLaunchDate,'Value',1);
end

% --- Executes during object creation, after setting all properties.
function menuConfiguration_CreateFcn(~, ~, ~) %#ok<DEFNU>
end

% --- Executes on button press in chkMoonLabeling.
function chkMoonLabeling_Callback(~, ~, handles) %#ok<DEFNU>
    if strcmp(get(handles.btnGo,'String'), 'Go') || ...
       strcmp(get(handles.btnGo,'String'), 'Resume')
        drawBodies(handles);
        drawLabels(0,0);
    end
end

% --- Executes on selection change in menuLaunchDate.
function menuLaunchDate_Callback(hObject, ~, handles) %#ok<DEFNU>
    global spaceship;
    global bodies;
    global bodyData;
    bodies = [];
    
    contents = cellstr(get(hObject,'String'));
    date = contents{get(hObject,'Value')};
    for i = 1:size(bodyData,1)
        name = bodyData{i,1};
        name = name{1};
        mass = bodyData{i,2};
        rad =  bodyData{i,3};
        color = 'k';
        pos = readtable(strcat('body_data/',name,'.txt'));
        
        for j=1:length(pos.epoch)
            if strcmp(pos.epoch{j},date)
                x = [pos.x(j); pos.y(j); pos.z(j)];
                v = [pos.vx(j); pos.vy(j); pos.vz(j)];
                %convert to AU/year
                v = v*365.242;                
                % pName,pPos,pVel,pMass,pRadius,pColor
                bodies = [bodies Body(name, x, v, mass, rad, color)];
                break;
            end
        end
        %pos contains all the position data for the launch date
    end
    %cassini is -82 in horizons
    cassini_pos = [-4.093868409910618E-02;  9.227191469266025E-01;  1.883044443995243E-02];
    cassini_vel = [-1.659025634165134E-02; -3.745936603666906E-03;  5.813019362138117E-05];
    cassini_vel = cassini_vel * 365.242;  
    % radius of ~10m, if it's that close then it's screwed anyway
    spaceship = Body('Ship',cassini_pos,cassini_vel,0.000001,7e-11,'b');
    bodies = [bodies spaceship];
    global venus;
    for b = bodies
        if strcmp(b.Name,'Venus')
            venus = b;
        end
    end
    
    comVel = [0;0;0];
    for b = bodies
        if ~ strcmp(b.Name, 'Sun')
            comVel = comVel + b.vel * b.Mass;
        end
    end
    for b = bodies
       if strcmp(b.Name, 'Sun')
           b.vel = -1 * comVel / b.Mass;
       end
    end
    drawBodies(handles);
    drawLabels(0,0);

end

function menuLaunchDate_CreateFcn(hObject, eventdata, handles)
end
