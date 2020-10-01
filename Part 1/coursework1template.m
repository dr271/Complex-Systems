% MCMCS Coursework 1
% Luc Berthouze 2018-10-23
% Leave the below unchanged, 
clear all


%--------------Choose landscape--------------
simpleLandscape = true; 
%False for complex landscape

%----------Plots the chosen landscape---------
if simpleLandscape
    lowerBound = -2;
    upperBound = 2;
    ezmesh(@SimpleLandscape,[-2 2],[-2 2])
else %If using Complex
    lowerBound = -3;
    upperBound = 7;
    ezmesh(@ComplexLandscape,[-3 7],[-3 7])
end
zlabel("Fitness");

%-------------Choose Parameters-------------
NumSteps=50;
LRate=0.1;
MaxMutate=1;
pcolorPlot = false; %Can't be used w/ individual startPt
individualStartPt = false;


if individualStartPt
    %-------------------------Select Start Point for individual algorithm Instance------------------
    StartPt = [-1,2];
    %StartPt = randomStartPt(lowerBound, upperBound);
    disp("Start co-ords: ");
    disp(StartPt);
    %---------------------------------------Select Algorithm--------------------------------------
    [maxReached, iterations, points] = GradAscent(StartPt,NumSteps,LRate, simpleLandscape);
    %[maxReached, iterations, points] = HillClimb(StartPt,NumSteps,MaxMutate, simpleLandscape);
    drawPath(points, simpleLandscape);

else
    %----------------------Alternative: Systematic grid of starting positions----------------------
        allPoints = [];
        i = 1;
        gridDensity = 1;
        %1 - Starting point on each integer co-ordinate
        %2 - Starting point on every 0.5 co-oridnate
        for x = gridDensity*lowerBound:gridDensity*upperBound
            for y = gridDensity*lowerBound:gridDensity*upperBound
                StartPt = [x/gridDensity y/gridDensity];
                %-----------------------------Select Algorithm-------------------------------------
                %[maxReached, iterations, points] = GradAscent(StartPt,NumSteps,LRate, simpleLandscape);
                [maxReached, iterations, points] = HillClimb(StartPt,NumSteps,MaxMutate, simpleLandscape);
                if ~ pcolorPlot drawPath(points, simpleLandscape), end 
                allPoints = cat(1, allPoints, points);%Used for pcolor plot
                performance(i, 1) = maxReached;
                performance(i, 2) = iterations;
                i = i + 1;
            end
        end
end





%-----------------------Code to generate pcolor plots----------------------------
if pcolorPlot && ~individualStartPt
    if simpleLandscape plotSize = 40, else plotSize = 110, end
        colourMap = zeros(plotSize);
        %Map tally of points on to coordinate grid for colourmap
        for i = 1:size(allPoints,1)
            %round and scale coordinates appropriately
            x = round((allPoints(i,1)+ abs(lowerBound) +.1)*10);
            if x > plotSize x=plotSize; end
            y = round((allPoints(i,2)+ abs(lowerBound) +.1)*10);
            if y > plotSize y=plotSize; end
            colourMap(x,y) = colourMap(x,y) + 1;
        end
        %pcolor omits last row and col, add one of each to negate
        newcol = zeros(plotSize,1);
        newrow = zeros(1,plotSize+1);
        colourMap = [colourMap newcol];
        colourMap = [colourMap; newrow];
        c = pcolor(linspace(lowerBound,upperBound,plotSize+1),linspace(lowerBound,upperBound,plotSize+1),colourMap');
        c.FaceColor = 'interp';
        colorbar;
end


%-----------------------Display Appropriate Output------------------------------
if individualStartPt && simpleLandscape
    disp("Global Optimum Reached: " + maxReached);
    disp("Iterations: " + iterations);

elseif ~individualStartPt && simpleLandscape
    [successProbability, avgIterations] = evaluateSuccess(performance);
    disp("Success Probability: " + successProbability);
    disp("Avg iterations to reach global optima: " +  avgIterations);
    
elseif individualStartPt && ~simpleLandscape
    disp("Maximum height reached: " + maxReached);
    
elseif ~individualStartPt && ~simpleLandscape
    
    [successProbability, avgIterations] = evaluateSuccess(performance);
    disp("Success Probability: " + successProbability);
    disp("Avg iterations to reach global optima: " +  avgIterations);
    
%     disp("Maximum height reached: " + max(performance(:,1)));
%     disp((max(performance(:,1))/12.235) + " of global-minima");
end



%========================================================================================================================

%A function that draws a line showing the path of chosen algorithm
%Params: coords - 2d array of x,y coords of all points on a path
%        simpleLandscape - Boolean indicating landscape in use
function drawPath(coords, simpleLandscape)
    for x = 1:size(coords,1)
        coords(x,3) = calculateHeight(coords(x,1),coords(x,2), simpleLandscape);
    end
    hold on
    plot3(coords(:,1), coords(:,2), coords(:,3)+0.1, 'r');
end


%Params: optima - 2d array mapping boolean success to iterations needed
%Returns: - Probability of reaching global optima
%         - Average iterations required to do so
function [successProbability, avgIterations] = evaluateSuccess(optima)
    successes = 0;
    iterations = 0;
    for x = 1:size(optima,1)
        if optima(x,1) == 1
            successes = successes + 1;
            iterations = iterations + optima(x,2);
        end
    end
    successProbability = successes/size(optima,1);
    avgIterations = iterations/successes;
end


%Returns: startPt - Random start point within bounds 
function [startPt] = randomStartPt(lowerBound, upperBound) 
    startPt = lowerBound + rand(1,2)*(upperBound-lowerBound);
end


%Function to calculate height/fitness
%Params: x, y - coordinates
%        simpleLandscape - boolean indicating landscape being used
%Returns: Height or fitness at given point
function [height] = calculateHeight(x, y, simpleLandscape)
    if simpleLandscape 
        height = SimpleLandscape(x,y);
    else
        height = ComplexLandscape(x,y);
    end
end


% Function implementing gradient ascent
%Returns: MaxReached - On simple landscape, binary if global optima found
%                    - On complex landscape, max height reached as double
%         Iterations - For simple landscape, iterations to reach G optima
function [maxReached, iterations, points] = GradAscent(StartPt,NumSteps,LRate, simpleLandscape)
    points = [];
	PauseFlag=0;
	hold on;
    maxReached = 0;
    
    
	for i = 1:NumSteps
        
	    %---------Calculates the 'height' at StartPt-----------------
        points = cat(1, points, StartPt);
        height = calculateHeight(StartPt(1),StartPt(2), simpleLandscape);
        disp("Height: " + height);
        %Updates max height for complex landscape output
        if height > maxReached
            maxReached = height;
%         elseif height == maxReached
%             disp("Stuck at local Maxima");
%             iterations = i;
%             break
        end
        
	    %-----------------Plots point on the landscape-----------------
        plot3(StartPt(1), StartPt(2), height, '*r', 'MarkerSize',10);
        
	
	    %------------Calculates the gradient at StartPt-----------------
        if simpleLandscape
            grad = SimpleLandscapeGrad(StartPt(1),StartPt(2));
        else
            grad = ComplexLandscapeGrad(StartPt(1),StartPt(2));
        end
        %disp("Gradient vector at current pos: ");
        %disp(grad);
        
	    %--------Calculates the new point and update StartPt-----------
        StartPt(1) = StartPt(1) + LRate*grad(1);
        StartPt(2) = StartPt(2) + LRate*grad(2);
        
            
	    %---------Ensures StartPt is within the specified bounds--------
        if simpleLandscape
            StartPt = max([StartPt;-2 -2]);
            StartPt = min([StartPt;2 2]);
        else
            StartPt = max([StartPt;-3 -3]);
            StartPt = min([StartPt;7 7]);
        end
            
        %--------------Checks if global optimum reached----------------
        %Overwrites maxReached to binary val if simple landscape
        iterations = i;
        if height == 4 && simpleLandscape || height >= 11.593705 && ~simpleLandscape
            disp("Global Optimum Reached")
            maxReached = 1;
            break
            %If step limit hit, return boolean 0 for maxima not met
        elseif i == NumSteps %%&& simpleLandscape
            maxReached = 0;
            break
%         elseif i == NumSteps && ~simpleLandscape
%             maxReached = i;
        end
        
        
	    %-------------------Pauses to view output----------------------
% 	    if(PauseFlag)
% 	        j=input('Press return to continue\nor 0 and return to stop pausing\n');
% 	        if(j==0) PauseFlag=0; end;
%         end
    end
    
    
	hold off
end

% Mutation function
% Params: OldPt - x,y location vector describing current position
%         MaxMutate - the maximum amount by which the current pos can be
%         mutated
% Returns: NewPt- The new, mutated potential location point 
function[NewPt] = Mutate(OldPt,MaxMutate)
	%Selects a random MutDist in the range(-MaxMutate,MaxMutate)
    MutDist = (2*MaxMutate)*rand(1,1) - MaxMutate;
	%---------Mutates random element of OldPt by MutDist-----------
    randElement = randi([1 2],1,1);
    OldPt(randElement) = OldPt(randElement)+ MutDist;
    NewPt = OldPt;
end	

%Returns: MaxReached - On simple landscape, binary if global optima found
%                    - On complex landscape, max height reached as double
%         Iterations - For simple landscape, iterations to reach G optima
function [maxReached, iterations, points] =  HillClimb(StartPt,NumSteps,MaxMutate, simpleLandscape)
    points = [];
	PauseFlag=0;
	hold on;
    maxReached = 0;
    noImprovementCount =0;
    
	for i = 1:NumSteps
        
        %---------Calculates the 'height' at StartPt-----------------
        points = cat(1, points, StartPt);%Comment out for improved pcolor plot
        height = calculateHeight(StartPt(1),StartPt(2), simpleLandscape);
        disp("Height: " + height);
        %Updates max height for complex landscape output
        if height > maxReached
            maxReached = height;
        end
        
        %---------------Plots point on landscape---------------------
        plot3(StartPt(1), StartPt(2), height, '*r', 'MarkerSize',10);
        
	    %-------------Mutates StartPt in to NewPt-------------------
	    NewPt=Mutate(StartPt,MaxMutate);

        %------Ensures NewPt is within the specified bounds---------
	    if simpleLandscape
            NewPt = max([NewPt;-2 -2]);
            NewPt = min([NewPt;2 2]);
        else
            NewPt = max([NewPt;-3 -3]);
            NewPt = min([NewPt;7 7]);
        end	

        %------Calculates the height of the new point--------------
        newHeight = calculateHeight(NewPt(1),NewPt(2), simpleLandscape);
        disp("New Height: " + newHeight);
                
        %-------Decides whether to update StartPt or not-----------   
        if newHeight >= height
            %points = cat(1, points, StartPt); %un-comment for improved pcolor plot
            StartPt = NewPt;
            disp("Point updated");
        else
            noImprovementCount = noImprovementCount + 1;
        end
        
        %----------------Terminating Conditions--------------------
        %Overwrites maxReached to binary val if simple landscape
        iterations = i;
        if height == 4 && simpleLandscape || height >= 11.593705 && ~simpleLandscape %Checks for global optimum
            disp("Global Optimum Reached")
            maxReached = 1;
            break
        elseif i == NumSteps %%&& simpleLandscape %Checks if iteration limit hit
            maxReached = 0;
            break
        elseif noImprovementCount == 50000
            maxReached = 0;
            disp("Stuck in local Minima");
            break
        end
        
	    %-----------------Pauses to view output-----------------
% 	    if(PauseFlag)
% 	        x=input('Press return to continue\nor 0 and return to stop pausing\n');
% 	        if(x==0) PauseFlag=0; end;
% 	    end		
	end
	hold off
end

% Definition of Simple landscape
function [z] = SimpleLandscape(x,y)
	z=max(1-abs(2*x),0)+x+y;
end

% Definition of gradient of Simple landscape
function [g] = SimpleLandscapeGrad(x,y) %Global Optima: 4
	if(1-abs(2*x) > 0)
	    if(x<0) g(1) = 3;
	    elseif(x==0) g(1)=0;
	    else g(1) = -1;
	    end
	else g(1)=1;
	end
	g(2)=1;
end

% Definition of Complex landscape
function [f]=ComplexLandscape(x,y) %Global Optima: 12.2039
	f=4*(1-x)^2*exp(-(x^2)-(y+1)^2) -15*(x/5 - x^3 - y^5)*exp(-x^2-y^2) -(1/3)*exp(-(x+1)^2 - y^2)-1*(2*(x-3)^7-0.3*(y-4)^5+(y-3)^9)*exp(-(x-3)^2-(y-3)^2);
end

% Definition of gradient of Complex landscape
function [g]=ComplexLandscapeGrad(x,y)
	g(1)=-8*exp(-(x^2)-(y+1)^2)*((1-x)+x*(1-x)^2)-15*exp(-x^2-y^2)*((0.2-3*x^2) -2*x*(x/5 - x^3 - y^5)) +(2/3)*(x+1)*exp(-(x+1)^2 - y^2)-1*exp(-(x-3)^2-(y-3)^2)*(14*(x-3)^6-2*(x-3)*(2*(x-3)^7-0.3*(y-4)^5+(y-3)^9));
	g(2)=-8*(y+1)*(1-x)^2*exp(-(x^2)-(y+1)^2) -15*exp(-x^2-y^2)*(-5*y^4 -2*y*(x/5 - x^3 - y^5)) +(2/3)*y*exp(-(x+1)^2 - y^2)-1*exp(-(x-3)^2-(y-3)^2)*((-1.5*(y-4)^4+9*(y-3)^8)-2*(y-3)*(2*(x-3)^7-0.3*(y-4)^5+(y-3)^9));
end
