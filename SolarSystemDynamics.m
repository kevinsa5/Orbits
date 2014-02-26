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

    % Last Modified by GUIDE v2.5 26-Feb-2014 16:37:44

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
    labeling = get(handles.chkLabeling, 'Value');
    lims = axis;
    cla;
    plot(spaceship.xHist,spaceship.yHist);
    hold on;
    for i=1:length(bodies)
        % using rectangle allows doing circle diameter in AU instead of px
        rectangle('Position',[transpose(bodies(i).pos - bodies(i).Radius) 2*bodies(i).Radius 2*bodies(i).Radius],'Curvature',1);
        plot(bodies(i).xHist, bodies(i).yHist, bodies(i).Color);
    end
    axis(lims);
    if labeling
        yax = ylim;
        xax = xlim;
        dx = 0.02*(xax(2)-xax(1));
        dy = 0.02*(yax(2)-yax(1));
        for i=1:length(bodies)
            %fprintf(1,'dx: %f; dy: %f \n', dx, dy);
            text(bodies(i).pos(1)+dx, bodies(i).pos(2)+dy, bodies(i).Name);
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
                    if i == j
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
                kr(:,1,i) = bodies(i).vel;
                accel = 0;
                for j = 1:length(bodies)
                    if i ~= j 
                        accel = accel + forceOn(bodies(i),bodies(j));
                    end
                end
                kv(:,1,i) = accel;
            end
            % calculate kr2, kv2 for each
            for i = 1:length(bodies)
                kr(:,2,i) = bodies(i).vel + kv(:,1,i)*deltaT/2;
                accel = 0;
                for j = 1:length(bodies)
                    if i ~= j
                        accel = accel + forceOn( bodies(i).rkCopy(deltaT/2 * kr(:,1,i)), ...
                                                 bodies(j).rkCopy(deltaT/2 * kr(:,1,j)));
                    end
                end
                kv(:,2,i) = accel;
            end
            % calculate kr3, kv3 for each 
            for i = 1:length(bodies)
                kr(:,3,i) = bodies(i).vel + kv(:,2,i)*deltaT/2;
                accel = 0;
                for j = 1:length(bodies)
                    if i ~= j
                        accel = accel + forceOn( bodies(i).rkCopy(deltaT/2 * kr(:,2,i)), ...
                                                 bodies(j).rkCopy(deltaT/2 * kr(:,2,j)));
                    end
                end
                kv(:,3,i) = accel;
            end
            % calculate kr4, kv4 for each 
            for i = 1:length(bodies)
                kr(:,4,i) = bodies(i).vel + kv(:,3,i)*deltaT;
                accel = 0;
                for j = 1:length(bodies)
                    if i ~= j
                        accel = accel + forceOn( bodies(i).rkCopy(deltaT * kr(:,3,i)), ...
                                                 bodies(j).rkCopy(deltaT * kr(:,3,j)));
                    end
                end
                kv(:,4,i) = accel;
            end
            for i=1:length(bodies)
                bodies(i).vel = bodies(i).vel + (deltaT/6) * (kv(:,1,i) + 2*kv(:,2,i) + 2*kv(:,3,i) + kv(:,4,i));
                bodies(i).pos = bodies(i).pos + (deltaT/6) * (kr(:,1,i) + 2*kr(:,2,i) + 2*kr(:,3,i) + kr(:,4,i));
            end
            for i=1:length(bodies)
                bodies(i).xHist = [bodies(i).xHist bodies(i).pos(1)];
                bodies(i).yHist = [bodies(i).yHist bodies(i).pos(2)];
            end
        elseif strcmp(method, 'Verlet')
            % Taken from http://xbeams.chem.yale.edu/~batista/vaa/node60.html
            dX = zeros(2,length(bodies));
            dV = zeros(2,length(bodies));
            for i = 1:length(bodies)
                accel = [0;0];
                iAccel = [0;0];
                for j = 1:length(bodies)
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
                bodies(i).pos = bodies(i).pos + dX(:,i);
                bodies(i).vel = bodies(i).vel + dV(:,i);

                bodies(i).xHist = [bodies(i).xHist bodies(i).pos(1)];
                bodies(i).yHist = [bodies(i).yHist bodies(i).pos(2)];
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
    % 10m = 6.685e-11 AU
    spaceship = Body331('Ship',4,0,0.001,6.685e-11,'k');
    sun = Body331('Sun', 0, 0, 333000, 0.004649, 'y');
    mercury = Body331('Mercury', 0.387,10.1,0.055,1.6308e-5,'r');
    venus = Body331('Venus', 0.723, 7.378, 0.8150, 4.0454e-5, 'g');
    earth = Body331('Earth', 1, 6.282, 1, 4.2564e-5, 'b');
    mars = Body331('Mars', 1.524, 5.076, 0.107, 2.263e-5,'r');
    jupiter = Body331('Jupiter', 5.203, 2.762, 317.8, 4.6239e-4, 'y');
    saturn = Body331('Saturn', 9.529, 2.02, 95.3, 3.8313e-4, 'y');
    uranus = Body331('Uranus', 19.19, 1.43, 14.6, 1.6889e-4, 'g');
    neptune = Body331('Neptune', 30.06, 1.14, 17.23, 1.6412e-4, 'b');

    bodies = [];
    
    if strcmp(configuration, 'Full Solar System')
        bodies = [spaceship sun mercury venus earth mars jupiter saturn uranus neptune];    
    elseif strcmp(configuration, 'Sun, Planets only')
        bodies = [spaceship sun mercury venus earth mars jupiter saturn uranus neptune];    
    elseif strcmp(configuration, 'Sun, Earth, Moon')
        bodies = [spaceship sun earth];
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

% --- Executes on selection change in menuMethod.
function menuMethod_Callback(~, ~, ~) %#ok<DEFNU>
    % Hints: contents = cellstr(get(hObject,'String')) returns menuMethod contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from menuMethod
end

% --- Executes during object creation, after setting all properties.
function menuMethod_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
    % Hint: popupmenu controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

% --- Executes on button press in chkLabeling.
function chkLabeling_Callback(~, ~, handles) %#ok<DEFNU>
% Hint: get(hObject,'Value') returns toggle state of chkLabeling
    if strcmp(get(handles.btnGo,'String'), 'Go') || ...
       strcmp(get(handles.btnGo,'String'), 'Resume')
        drawBodies(handles);
    end
end

function menuConfiguration_Callback(hObject, ~, handles) %#ok<DEFNU>
    global bodies;
    contents = cellstr(get(hObject,'String'));
    config = contents{get(hObject,'Value')};
    bodies = DefineBodies(config);
    drawBodies(handles);
end