function varargout = SolarSystemDynamics(varargin)
    % SOLARSYSTEMDYNAMICS MATLAB code for SolarSystemDynamics.fig
    %      SOLARSYSTEMDYNAMICS, by itself, creates a new SOLARSYSTEMDYNAMICS or raises the existing
    %      singleton*.
    %
    %      H = SOLARSYSTEMDYNAMICS returns the handle to a new SOLARSYSTEMDYNAMICS or the handle to
    %      the existing singleton*.
    %
    %      SOLARSYSTEMDYNAMICS('CALLBACK',hObject,~,handles,...) calls the local
    %      function named CALLBACK in SOLARSYSTEMDYNAMICS.M with the given input arguments.
    %
    %      SOLARSYSTEMDYNAMICS('Property','Value',...) creates a new SOLARSYSTEMDYNAMICS or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before SolarSystemDynamics_OpeningFcn gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to SolarSystemDynamics_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % Edit the above text to modify the response to help SolarSystemDynamics

    % Last Modified by GUIDE v2.5 29-Mar-2014 16:43:54

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
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % ~  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to SolarSystemDynamics (see VARARGIN)

    % Choose default command line output for SolarSystemDynamics
    handles.output = hObject;
    % Update handles structure
    guidata(hObject, handles);
    % UIWAIT makes SolarSystemDynamics wait for user response (see UIRESUME)
    % uiwait(handles.MainFigure);
    reset(handles);
    menuConfiguration_Callback(hObject,0,handles);
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


% --- Executes on slider movement.
function initial_velocity_Callback(hObject, ~, handles) %#ok<DEFNU>
    value = get(hObject,'Value');
    set(handles.txt_velocity,'String',value);
    initial_angle_Callback(handles.initial_angle, 0, handles);
end

function drawBodies(handles)
    global bodies
    global spaceship;
    planetLabeling = get(handles.chkLabeling, 'Value');
    moonLabeling = get(handles.chkMoonLabeling, 'Value');
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
        rectangle('Position',[transpose(bodies(i).pos - bodies(i).Radius) 2*bodies(i).Radius 2*bodies(i).Radius],'Curvature',1);
        plot(bodies(i).xHist, bodies(i).yHist, bodies(i).Color);
    end
    % used for centering on the currently tracked body:
    deltaAxis = [translation(1) - mean(lims(1:2)), translation(2) - mean(lims(3:4))];
    axis(lims + [deltaAxis(1) deltaAxis(1) deltaAxis(2) deltaAxis(2)]);
    
    if planetLabeling || moonLabeling
        yax = ylim;
        xax = xlim;
        dx = 0.02*(xax(2)-xax(1));
        dy = 0.02*(yax(2)-yax(1));
        for i=1:length(bodies)
            if bodies(i).joined, continue; end
            if (strcmp(bodies(i).Name,'Sun') || strcmp(bodies(i).parent.Name,'Sun')) && planetLabeling
                % bodies(i) is a planet or the sun
                text(bodies(i).pos(1)+dx, bodies(i).pos(2)+dy, bodies(i).Name);
            elseif bodies(i).numberOfChildren == 0 && ~ strcmp(bodies(i).parent.Name,'Sun') && moonLabeling
                % bodies(i) is a moon
                text(bodies(i).pos(1)+dx, bodies(i).pos(2)-dy*1.5*(bodies(i).childNumber), bodies(i).Name);
            end
            
        end
    end
    hold off;
end

function btnGo_Callback(hObject, ~, handles) %#ok<DEFNU>
    global bodies;
    global spaceship;

    if strcmp(get(hObject, 'String'), 'Pause')
        disp('pause was hit');
        set(hObject, 'String', 'Resume');
        set(handles.btnReset, 'Enable','on');
        return;
    end

    disp('started');
    start = str2double(get(hObject,'UserData'));
    if start == 1 
        v0 = get(handles.initial_velocity,'value');
        t0 = get(handles.initial_angle,'value');
        spaceship.vel = [v0*cos(t0); v0*sin(t0)];
    end
    set(hObject, 'String', 'Pause');
    set(handles.initial_velocity,'Enable','off');
    set(handles.initial_angle, 'Enable', 'off');
    set(handles.menuConfiguration, 'Enable', 'off');
    set(handles.txtTimeStep, 'Enable','off');
    set(handles.txtFrameSkip, 'Enable','off');
    set(handles.btnReset, 'Enable','off');
    set(handles.menuMethod, 'Enable', 'off');
    
    deltaT = str2double(get(handles.txtTimeStep, 'String'));
    frameSkip = str2double(get(handles.txtFrameSkip, 'String'));

    for i = bodies
        fprintf(1,'Name: %s; x: %f; y: %f;\n', i.Name, i.pos(1), i.pos(2));
    end

    % 10000 is the arbitrary number of iterations I picked
    for t=start:70000
        % window was closed?
        if ~ ishandle(handles.MainFigure)
            break;
        end

        if strcmp(get(hObject, 'String'), 'Resume')
            set(hObject, 'UserData', num2str(t));
            %drawBodies(true);
            drawBodies(handles);
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
            kr = zeros(2,order,length(bodies));
            kv = zeros(2,order,length(bodies));
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
            dX = zeros(2,length(bodies));
            dV = zeros(2,length(bodies));
            for i = 1:length(bodies)
                if bodies(i).joined, continue; end
                accel = [0;0];
                iAccel = [0;0];
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

        %axis(defaultAxes);
        if mod(t,frameSkip) == 0 
            drawBodies(handles);
            drawnow;
        end

    end
    if ishandle(hObject)
        set(hObject, 'UserData', 'Finished');
        set(hObject, 'String', 'Go');
        set(hObject, 'UserData', '1');
    end
end

function accel = forceOn(body1, body2)
    % Big G in units of Au^3 / (earth mass * year^2)
    G = 0.00011835;
    r = norm(body1.pos - body2.pos);
    % next line has (x1-x2)/mag(x1-x2) to represent rhat
    force = -G*body1.Mass*body2.Mass/ r^3 * (body1.pos - body2.pos);
    accel = force / body1.Mass;
end
    
function initial_angle_Callback(hObject, ~, handles)
    global spaceship;
    value = get(hObject,'Value');
    set(handles.txt_angle,'String',value*180/pi);
    if ~ strcmp(get(handles.btnGo, 'String'), 'Go')
        return;
    end
    tempr = 0:0.01:get(handles.initial_velocity, 'Value');
    drawBodies(handles);
    hold on;
    plot(spaceship.pos(1)+tempr*cos(value),spaceship.pos(2)+tempr*sin(value));
    %axis(defaultAxes);
    hold off;
end

function bodies = DefineBodies(configuration)
    global spaceship;
    global sun;
    % pName,pOrbitRadius,pOrbitSpeed,pMass,pRadius,pColor
    % 10m = 6.685e-11 AU
    sun = Body('Sun', 0, 0, 333000, 0.004649, 'y');
        spaceship = sun.makeMoon('Ship',3.8,2,0.001,6.685e-11,'k');
        mercury = sun.makeMoon('Mercury', 0.387,10.1,0.055,1.6308e-5,'r');
        venus = sun.makeMoon('Venus', 0.723, 7.378, 0.8150, 4.0454e-5, 'g');
        earth = sun.makeMoon('Earth', 1, 6.282, 1, 4.2564e-5, 'b');
            earthMoon = earth.makeMoon('Mona', 0.00257, 0.21544, 0.012300, 1.161e-5,'k');
        mars = sun.makeMoon('Mars', 1.524, 5.076, 0.107, 2.263e-5,'r');
            phobos = mars.makeMoon('Phobos', 6.2675e-5, 0.4507, 1.78477e-8, 7.5e-8,'k');
            deimos = mars.makeMoon('Deimos',1.5684E-004, 2.8490E-002, 2.4718E-008, 4.1000E-008,'k');
        jupiter = sun.makeMoon('Jupiter', 5.203, 2.762, 317.8, 4.6239e-4, 'y');
            io = jupiter.makeMoon('Io',2.8190E-003, 3.6540E+000, 1.5000E-002, 1.2175E-005,'k');
            europa = jupiter.makeMoon('Europa',4.4847E-003, 2.8960E+000, 8.0000E-003, 1.0433E-005,'k');
            ganymede = jupiter.makeMoon('Ganymede', 7.1552E-003, 2.2940E+000, 2.5000E-002, 1.7608E-005,'k');
            callisto = jupiter.makeMoon('Callisto', 1.2585E-002, 1.7294E+000, 1.8000E-002, 1.6110E-005,'k');
        saturn = sun.makeMoon('Saturn', 9.529, 2.02, 95.3, 3.8313e-4, 'y');
            mimas = saturn.makeMoon('Mimas', 1.2366E-003, 3.0107E+000, 6.3000E-006, 1.3250E-006,'k');
            enceladus = saturn.makeMoon('Enceladus', 1.5906E-003, 2.6622E+000, 1.8000E-005, 1.6850E-006,'k');
            tethys = saturn.makeMoon('Tethys', 1.9694E-003, 2.3925E+000, 1.0300E-004, 3.5500E-006,'k');
            dione = saturn.makeMoon('Dione', 2.5227E-003, 2.1141E+000, 3.2800E-004, 3.7530E-006,'k');
            rhea = saturn.makeMoon('Rhea', 3.5235E-003, 1.7885E+000, 3.9000E-004, 5.1060E-006,'k');
            titan = saturn.makeMoon('Titan', 8.1677E-003, 1.1748E+000, 2.2500E-002, 1.7219E-005,'k');
            iapetus = saturn.makeMoon('Iapetus', 2.3803E-002, 6.8820E-001, 3.0234E-004, 4.9100E-006,'k');
        uranus = sun.makeMoon('Uranus', 19.19, 1.43, 14.6, 1.6889e-4, 'g');
            miranda = uranus.makeMoon('Miranda', 8.6492E-004, 0.223346, 1.1030E-005, 1.5760E-006,'k');
            ariel = uranus.makeMoon('Ariel', 1.2769E-003, 1.8495E-001, 2.2600E-004, 3.8700E-006,'k');
            umbriel = uranus.makeMoon('Umbriel', 1.7781E-003, 0.156614, 2.0000E-004, 3.9080E-006,'k');
            tatania = uranus.makeMoon('Titania', 2.9139E-003, 1.2216E-001, 5.9080E-004, 5.2700E-006,'k');
            oberon = uranus.makeMoon('Oberon', 3.9006E-003, 0.105748, 5.0460E-004, 5.0900E-006,'k');
        neptune = sun.makeMoon('Neptune', 30.06, 1.14, 17.23, 1.6412e-4, 'b');
            triton = neptune.makeMoon('Triton', 2.3714E-003, 0.147284, 3.5900E-003, 9.0470E-006,'k');
        
    bodies = [];
    
    if strcmp(configuration, 'Full Solar System')
        bodies = [spaceship sun mercury venus earth mars jupiter saturn uranus neptune earthMoon phobos deimos io europa ganymede callisto mimas enceladus tethys dione rhea titan iapetus miranda ariel umbriel tatania oberon triton];    
    elseif strcmp(configuration, 'Sun, Planets only')
        bodies = [spaceship sun mercury venus earth mars jupiter saturn uranus neptune];    
    elseif strcmp(configuration, 'Sun, Earth, Moon')
        bodies = [spaceship sun earth earthMoon];
    end
    
    % make the system's net momentum zero:
    comVel = [0;0];
    for b = bodies
        if ~ strcmp(b.Name, 'Sun')
            comVel = comVel + b.vel * b.Mass;
        end
    end
    sun.vel = -1 * comVel / sun.Mass;
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
    global bodies;

    defaultAxes = [-30 30 -30 30];
    axis(defaultAxes);
    set(handles.txtTimeStep, 'String', 0.00273);
    set(handles.txtFrameSkip, 'String', 10);
    set(handles.initial_angle, 'Value', pi/2);
    set(handles.initial_velocity,'Value',1);
    set(handles.btnGo, 'UserData', '1');
    set(handles.txt_velocity,'String',get(handles.initial_velocity,'Value'));
    set(handles.txt_angle,'String',180/pi*get(handles.initial_angle,'Value'));
    set(handles.FrameCount,'String',0);
    
    set(handles.initial_velocity,'Enable','on');
    set(handles.initial_angle, 'Enable', 'on');
    set(handles.menuConfiguration, 'Enable', 'on');
    set(handles.txtTimeStep, 'Enable','on');
    set(handles.txtFrameSkip, 'Enable','on');
    set(handles.menuMethod, 'Enable', 'on');
    set(handles.btnGo, 'String', 'Go');
    
    contents = cellstr(get(handles.menuConfiguration,'String'));
    config = contents{get(handles.menuConfiguration,'Value')};
    bodies = DefineBodies(config);
    drawBodies(handles);
end

function btnReset_Callback(~, ~, handles) %#ok<DEFNU>
    reset(handles);
end

function menuMethod_Callback(~, ~, handles) %#ok<DEFNU>
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
    end
end

function menuConfiguration_Callback(~, ~, handles) %#ok<DEFNU>
    global bodies;
    contents = cellstr(get(handles.menuConfiguration,'String'));
    config = contents{get(handles.menuConfiguration,'Value')};
    bodies = DefineBodies(config);
    nameList{1} = 'None';
    for i = 1:length(bodies)
        nameList{i+1} = bodies(i).Name;
    end
    set(handles.centeredBody,'String',nameList);
    set(handles.centeredBody,'Value',1);
    drawBodies(handles);   
end


% --- Executes during object creation, after setting all properties.
function menuConfiguration_CreateFcn(~, ~, ~) %#ok<DEFNU>
end


% --- Executes on button press in chkMoonLabeling.
function chkMoonLabeling_Callback(~, ~, handles) %#ok<DEFNU>
    if strcmp(get(handles.btnGo,'String'), 'Go') || ...
       strcmp(get(handles.btnGo,'String'), 'Resume')
        drawBodies(handles);
    end
end