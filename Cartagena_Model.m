clear all; close all; clc;
% Number of iterationa
numIter = 1500;
transitorio=500;

% Values of h to simulate
h = 0.8;

% Load Network information
[numDataNode,txtDataNode] = xlsread('Nodos_Ctg_1.xlsx',1);
[numDataLink,txtDataLink] = xlsread('vias_copia_1.xlsx',1);

% Size of the network
N = length(numDataNode);

% XY coordinates of the nodes
spaceXYnode = numDataNode(:,1:2);
nodeId = numDataNode(:,3);

% Number of links
L = length(numDataLink);

% Initial and Final Node of the links
nodoIni = numDataLink(:,8);
nodoFin = numDataLink(:,9);
% Length of the link
LengthLink = numDataLink(:,16);

% Distance matrix
D = zeros(N);
for k = 1:L
    i = find(nodeId == nodoIni(k));
    j = find(nodeId == nodoFin(k));
    D(i,j) = LengthLink(k);
    D(j,i) = LengthLink(k);
end

% Eliminate links with only 1 intersection
[Dcoarse,spCoarseXY] = coarseGrainGraph(D,spaceXYnode);

% Size of the Coarse Grained graph
N = length(Dcoarse);

A = zeros(N);
A(Dcoarse>0)=sparse(1./Dcoarse(Dcoarse>0));

% Construct graph using distances as weights of the link
G = graph(Dcoarse);

% Shortest-path Matrix
D = distances(G);
D = D/mean(Dcoarse(Dcoarse>0));

% R parameter
R = 5;

% Transfer capacity per step
C = 1;

% Maximum capacity of the network:
maxPaquetes = R*(numIter+1);

% Vector with agent Identification
agId = zeros(maxPaquetes,1);

% Matrix with information of initial node, destination node and current node of each agent
agState = zeros(maxPaquetes,3);

% Waitibg time histogram parameters
maxWaiting = 300; minWaiting = 1; nBins = 50;
waitingTime = zeros(maxPaquetes,1);

waitingBins = linspace(minWaiting,maxWaiting,nBins);
waitingHist = zeros(1,nBins);

% Cell storing the agents in the queue of the node
cuePacket = cell(N,1);

% Initilize the number of agents in the network
Np_old = R;

% Initialize order parameter
order_parameter = zeros(numIter,1);
numbOfPack = zeros(numIter,1);

% Initial occupation
ocup = zeros(N,1);
IdxVec = 1:maxPaquetes;      % At the beginning all Identifiers are Free

freeIDX = find(IdxVec>0);    % Identify free Identifiers
auxIDX = freeIDX(1:R);       % Take the first R free Idnetifiers
agId(auxIDX) = auxIDX;       % Assign the Id of the new agents in the network the free Identifiers
IdxVec(auxIDX)=0;            % The assigned identifiers are not free anymore, set to 0

ini_node = randi(N,[R 1]);   % Create randomly R initial nodes
end_node = randi(N,[R 1]);   % Create randomly R destination nodes
agState(auxIDX,1:2) = [ini_node end_node]; % Assign initial node and destination node to first columns of agState matrix
agState(auxIDX,3) = ini_node; % Current node is the initial one

for j = 1:N
    cuePacket{j} = find(agState(:,3)==j); % Update queue finding which of the agents in the network are in each node
end


% Draw the layout of the map from the shape files

S = shaperead('SHAPES/VIAS.shp');

figure(4)
hold on
F = mapshow(S);
set(gca,'Visible','on','YLim',[1.635e6 1.65e6],'XLim',[8.36e5 8.52e5])

h1 = scatter(spCoarseXY(:,1),spCoarseXY(:,2),20,'filled');
colormap(jet(20))
set(gca,'Clim',[0 20]);
colorbar

% Iterate the model
for i = 1:numIter            
    
    if(mod(i,100)==0)
        i
    end
    
    for j = 1:N
        % For each node extract the queue and the neighbors of that node
        cue = cuePacket{j};
        neighbors = find(A(j,:)>0);
        
        nIni = j;
        
        itNode = min(C,length(cue)); % Number of agents to transfer, in the paper this is 1.
        for k = 1:itNode
            id = cue(k); % Extract the first agent of the queue
            nFin = agState(id,2); % Check the destination of that agent
            % Check if one of the neighbors to the node is destination 
            if(~ismember(nFin,neighbors)) % If not:
                neigDistance = D(neighbors,nFin);  % Calculate distance between the possible nodes and destination node
                neigCue = zeros(size(neigDistance)); % Initialize vector where I will store the size of queue of the neighbors
                
                for l = 1:length(neighbors)
                    neigCue(l) = length(cuePacket{neighbors(l)}); % Size of the queue of the neighbors to the node
                end
                
                
                delta_i = h*neigDistance + (1-h)*neigCue;  % effective distance
                
                idNeigh = find(delta_i==min(delta_i)); % Find the neighbor minimizing the effective distance
                
                if(length(idNeigh)==1) % If only one neighbor meets the cirteria of minimal distance, move the agent to that node
                    agState(id,3)=neighbors(idNeigh);
                else
                    idNeigh = randsample(idNeigh,1); % If there is more than one choose randmoly between the neighbors with minimum effective distance
                    agState(id,3)=neighbors(idNeigh); 
                end
                
                waitingTime(id)=waitingTime(id)+1; % we add 1 iteration to the waiting time of that agent
                
                % If one of the neighbors was the destination node, then:
            else
                agState(id,3)=-1; % The agent arrived to destination, state is set to -1
                IdxVec(id) = id; % The arrived agent sets free its identifier that was 0 during the time it was moving through the network
                
                % This updates the information of the waiting time
                % histogram:
                if(waitingTime(id)~=0) 
                    binWaiting = floor(waitingTime(id)*nBins/(maxWaiting-minWaiting))+1; 
                    if (binWaiting<=nBins && i>transitorio)
                        waitingHist(binWaiting)=waitingHist(binWaiting)+1; 
                    end
                    waitingTime(id)=0;
                    
                end
                
            end
            
        end
        
    end
    
    Np_new = sum(agState(:,3)>0); % Number of agents in current iteration, are those with destination state larger than 0
    deltaNp = Np_new-Np_old; % Difference between current number of agents and previous iteration
    numbOfPack(i) = Np_new; % Store the number of agents in this iteration
    order_parameter(i) = deltaNp; % Order parameter
    Np_old = Np_new; % Update current number of agents
    
    % As in the beginning this assings identifiers, origin and destination nodes
    % for the next R agents created
    freeIDX = find(IdxVec>0);
    auxIDX = freeIDX(1:R);
    agId(auxIDX) = auxIDX;
    IdxVec(auxIDX)=0;
    
    ini_node = randi(N,[R 1]);
    end_node = randi(N,[R 1]);
    agState(auxIDX,1:2) = [ini_node end_node];
    agState(auxIDX,3) = ini_node;
    
    for j = 1:N
        cuePacket{j} = find(agState(:,3)==j);
        ocup(j)=length(cuePacket{j});
    end
    
    % Update graphic showing number of agents in each node
    h1.CData = ocup;
    drawnow        
          
end

%% Draw order parameter and occupation of the network

figure
subplot(2,1,1)
plot(numbOfPack)
xlabel('iteration'); ylabel('load');
subplot(2,1,2)
plot(order_parameter)
xlabel('iteration'); ylabel('\eta');

%% Draws the distribution of waiting times
v=waitingHist/trapz(waitingBins,waitingHist);

figure
plot(waitingBins,v,'o')
set(gca,'XScale','lin','YScale','log')
xlabel('Tiempo de espera'); ylabel('Probabilidad')
parametrodeorden=(numbOfPack(end)-numbOfPack(transitorio))/((numIter-transitorio)*R);
