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

% Last Modified by GUIDE v2.5 17-Feb-2014 08:45:18

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

% --- Outputs from this function are returned to the command line.
function varargout = SolarSystemDynamics_OutputFcn(hObject, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function initial_velocity_Callback(hObject, ~, handles)
% hObject    handle to initial_velocity (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
value = get(hObject,'Value');
set(handles.txt_velocity,'String',value);
initial_angle_Callback(handles.initial_angle, 0, handles);

function drawBodies(varargin)
global bodies
global spaceship;
if nargin > 0
    labeling = varargin{1};
else
    labeling = false;
end
lims = axis;
cla;
plot(spaceship.xHist,spaceship.yHist);
hold on;
for i=1:length(bodies)
    scatter(bodies(i).pos(1), bodies(i).pos(2),bodies(i).Mass^(1/9)*100,strcat(bodies(i).Color,'*'));
    plot(bodies(i).xHist, bodies(i).yHist, bodies(i).Color);
    if labeling
        yax = ylim;
        xax = xlim;
        dx = 0.05*(xax(2)-xax(1));
        dy = 0.05*(yax(2)-yax(1));
        %fprintf(1,'dx: %f; dy: %f \n', dx, dy);
        text(bodies(i).pos(1)+dx, bodies(i).pos(2)+dy, bodies(i).Name);
    end
end
axis(lims);
hold off;

function btnGo_Callback(hObject, ~, handles)
global bodies;
global spaceship;
global sun;
% hObject    handle to btnGo (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(hObject, 'String'), 'Pause')
    disp('pause was hit');
    set(hObject, 'String', 'Resume');
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
        drawBodies();
        %axis(defaultAxes);
        return;
    end
    set(hObject, 'UserData', strcat('Running ', num2str(t)));
    set(handles.FrameCount, 'String', num2str(t));

    
    if mod(t,frameSkip) == 0
        drawBodies();
    end
 
    contents = cellstr(get(handles.menuMethod,'String'));
    method = contents{get(handles.menuMethod,'Value')};
    if strcmp(method,'Forward Euler')
        %eventually, do j=(i+1):length(bodies) and accelerate j based on i as well.
        %as simple as just adding pi to theta?
        %right now, we're duplicating all Gm/r^2 calculations
        for i=1:length(bodies)
            for j=1:length(bodies)
                if i == j
                    continue
                end
                force = forceOn(bodies(i),bodies(j));
                a = force / bodies(i).Mass;
                bodies(i).vel = bodies(i).vel + a*deltaT;
            end
        end
        % now, move all the bodies
        % can't combine with the acceleration loop because you can't calculate
        % forces on one, then move it, then use that new position for other
        % bodies' forces. move all, then accelerate all.
        for i=1:length(bodies)
            bodies(i).pos = bodies(i).pos + bodies(i).vel * deltaT;
    %        if strcmp(bodies(i).Name, 'Sun')
    %            disp(bodies(i).pos(1)-sun.pos(1));
    %        end
            bodies(i).xHist = [bodies(i).xHist bodies(i).pos(1)];
            bodies(i).yHist = [bodies(i).yHist bodies(i).pos(2)];
        end
    elseif strcmp(method,'Runge Kutta 4')
        for i=1:length(bodies)
            for j=1:length(bodies)
                if i == j
                    continue
                end
                % see http://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#The_Runge.E2.80.93Kutta_method
                % This implementation is described in http://spiff.rit.edu/richmond/nbody/OrbitRungeKutta4.pdf
                kv = [[0;0] [0;0] [0;0] [0;0]];
                kr = [[0;0] [0;0] [0;0] [0;0]];
                
                
                kr(:,1) = bodies(i).vel;
                kv(:,1) = forceOn(bodies(i),bodies(j)) / bodies(i).Mass;
                kr(:,2) = bodies(i).vel .* kv(1) * deltaT/2;
                kv(:,2) = forceOn(bodies(i).rkCopy(kr(1)*deltaT/2),bodies(j)) / bodies(i).Mass;
                kr(:,3) = bodies(i).vel .* kv(2) * deltaT/2;
                kv(:,3) = forceOn(bodies(i).rkCopy(kr(2)*deltaT/2),bodies(j)) / bodies(i).Mass;
                kr(:,4) = bodies(i).vel .* kv(3) * deltaT;
                kv(:,4) = forceOn(bodies(i).rkCopy(kr(3)*deltaT),bodies(j)) / bodies(i).Mass;
                
                bodies(i).vel = bodies(i).vel + (deltaT/6) * (kv(1) + 2*kv(2) + 2*kv(3) + kv(4));
                bodies(i).pos = bodies(i).pos + (deltaT/6) * (kr(1) + 2*kr(2) + 2*kr(3) + kr(4));
            end
        end
        for i=1:length(bodies)
            bodies(i).xHist = [bodies(i).xHist bodies(i).pos(1)];
            bodies(i).yHist = [bodies(i).yHist bodies(i).pos(2)];
        end
    elseif strcmp(method, 'Verlet')
        a = zeros(2,length(bodies));
        for i=1:length(bodies)
            if size(bodies(i).xHist)==0
                hist = bodies(i).pos - bodies(i).vel * deltaT;
                bodies(i).xHist = hist(1);
                bodies(i).yHist = hist(2);
            end
            netForce = 0;
            for j=1:length(bodies)
                if i == j
                    continue
                end
                netForce = netForce + forceOn(bodies(i),bodies(j));
            end
            a(:,i) = netForce / bodies(i).Mass;
        end
        for i=1:length(bodies)
            bodies(i).pos = bodies(i).pos + bodies(i).vel * deltaT + 0.5 * a(i) * deltaT^2;
            bodies(i).vel = bodies(i).vel + 0.5 * (a(i) + ) * deltaT;
            %bodies(i).pos = bodies(i).pos + (bodies(i).pos - [bodies(i).xHist(end); bodies(i).yHist(end)]) + a(:,i) * deltaT * deltaT;
            bodies(i).xHist = [bodies(i).xHist bodies(i).pos(1)];
            bodies(i).yHist = [bodies(i).yHist bodies(i).pos(2)];
        end
    end

    %axis(defaultAxes);
    if mod(t,frameSkip) == 0 
        drawnow;
    end
    
end
if ishandle(hObject)
    set(hObject, 'UserData', 'Finished');
    set(hObject, 'String', 'Go');
    set(hObject, 'UserData', '1');
end


function force = forceOn(body1, body2)
    % Big G in units of Au^3 / (earth mass * year^2)
    %G = 3.964e29;
    G = 0.00011835;
    r = sqrt(sum((body1.pos-body2.pos).^2));
    f = -G*body1.Mass*body2.Mass/r^2;
    theta = atan((body1.pos(2)-body2.pos(2))/(body1.pos(1)-body2.pos(1)));
    %atan has a funny range, this turns it into 0-2pi
    if body1.pos(1) < body2.pos(1)
        theta = theta + pi;
    elseif body1.pos(2) < body2.pos(2)
        theta = theta + 2*pi;
    end
    force = [ f*cos(theta) ; f*sin(theta) ];
    
    
function initial_angle_Callback(hObject, ~, handles)
global spaceship;
value = get(hObject,'Value');
set(handles.txt_angle,'String',value*180/pi);
if ~ strcmp(get(handles.btnGo, 'String'), 'Go')
    return;
end
tempr = 0:0.01:get(handles.initial_velocity, 'Value');
drawBodies();
hold on;
plot(spaceship.x+tempr*cos(value),spaceship.y+tempr*sin(value));
%axis(defaultAxes);
hold off;

function bodies = DefineBodies()
global spaceship;
global sun;

configuration = 'Solar System';

if strcmp(configuration, 'Solar System')
    spaceship = Body331('Ship',50,0,0.001,0.1,'k');
    %  Body(   Name,  orbitRadius, orbitSpeed ,Mass,  Radius, Color)
    sun = Body331('Sun', 0, 0, 333000, 109, 'y');
    mercury = Body331('Mercury', 0.387,10.1,0.055,0.382 / 2,'r');
    venus = Body331('Venus', 0.723, 7.378, 0.8150, 0.949 / 2, 'g');
    earth = Body331('Earth', 1, 6.282, 1, 1, 'b');
    mars = Body331('Mars', 1.524, 5.076, 0.107, 0.532 / 2,'r');
    jupiter = Body331('Jupiter', 5.203, 2.762, 317.8, 11.19/2, 'y');
    saturn = Body331('Saturn', 9.529, 2.02, 95.3, 9.26/2, 'y');
    uranus = Body331('Uranus', 19.19, 1.43, 14.6, 4.01/2, 'g');
    neptune = Body331('Neptune', 30.06, 1.14, 17.23, 3.88 / 2, 'b');
   
    %earthMoon = Body('Earth''s Moon', 0, 50.1, -20, 0, 0.5, 4000, 'k');
    bodies = [spaceship sun mercury venus earth mars jupiter saturn uranus neptune];    
end

function txtTimeStep_Callback(hObject, ~, handles)
% hObject    handle to txtTimeStep (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtTimeStep as text
%        str2double(get(hObject,'String')) returns contents of txtTimeStep as a double
val = get(hObject, 'String');
if isnan(str2double(val))
    set(hObject, 'String', 0.01);
    msgbox(strcat('Cannot convert ''', val, ''' to a number'),'Time Step not valid');
elseif str2double(val) <= 0
    set(hObject, 'String', 0.01);
    msgbox('Time Step must be greater than zero', 'Time Step not valid');
end

function txtFrameSkip_Callback(hObject, ~, handles)
% hObject    handle to txtFrameSkip (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtFrameSkip as text
%        str2double(get(hObject,'String')) returns contents of txtFrameSkip as a double
val = get(hObject, 'String');
if isnan(str2double(val))
    set(hObject, 'String', 1);
    msgbox(strcat('Cannot convert ''', val, ''' to a number'), 'Frame Skip not valid');
elseif str2double(val) <= 0
    set(hObject, 'String', 1);
    msgbox('Frame Skip must be greater than zero', 'Frame Skip not valid');
end

function MainFigure_CloseRequestFcn(hObject, ~, handles)
% hObject    handle to MainFigure (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);

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

bodies = DefineBodies();
drawBodies();


function btnReset_Callback(hObject, ~, handles)
reset(handles);

% --- Executes on selection change in menuMethod.
function menuMethod_Callback(hObject, ~, handles)
% hObject    handle to menuMethod (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns menuMethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menuMethod


% --- Executes during object creation, after setting all properties.
function menuMethod_CreateFcn(hObject, ~, handles)
% hObject    handle to menuMethod (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
