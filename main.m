clear all;

load xy.mat;

iteration =1000;
population_size = 200;
optRoute = zeros(100,2);
population = zeros(population_size,100);
test=generatePath();
for i=1:population_size
    population(i,:)=generatePath();
    
end
%% add a new row to the population to be used for distance
population =[population zeros(population_size,1)];

%   P1=[4 7 1 6 5 2 3 8]
%   P2=[1 2 3 4 8 5 6 7]
%   
%   test = orderCross(P1,P2)

%%
for j=1:iteration
    %% calculate the distance for every individual
    for i=1:population_size
        population(i,101) = distance(population(i,:),xy);
    end
    
    %% elite
    population = sortrows(population,101);
    newPopulation = zeros(population_size,100);
    newPopulation(1:2,:)=population(1:2,1:100);
    copy=population(1:100,:);
    copy(:,end)=[];
    population_new_num=2;
    %% Create a new generation
    while(population_new_num<population_size+1)
        %% crossover 
        %always the best from the previous generation to a selected 
         parent1=population(1,1:100);
         parent1=select(copy,xy);
         parent2=select(copy,xy);
        %always cross over
        firstChild = crossOver(parent1,parent2);
        secondChild = crossOver(parent2,parent1);
        
        %% mutation
            if rand<0.5
                firstChild = IMutate(firstChild);
           end
           if rand<0.5
               secondChild = IMutate(secondChild);
           end
        %% store the offspring in a new generation
        %find the best and store it
        bestChild=best(firstChild,secondChild,xy);
        newPopulation(population_new_num,:)=bestChild;
        population_new_num = population_new_num+1;
    end
    population(:,1:100)=newPopulation;
end
route = population(1,:);
plotGraph(route,xy);
%% plot the data
%plot the best route
function plotGraph(optRoute,xy)
    figure('Name','TSP_GA | Results','Numbertitle','off');         
    subplot(2,2,1);         
    pclr = ~get(0,'DefaultAxesColor');         
    plot(xy(:,1),xy(:,2),'.','Color',pclr);          
    title('City Locations');         
    subplot(2,2,2);         
    rte = optRoute([1:100 1]);         
    plot(xy(rte,1),xy(rte,2),'r.-');         
    title(sprintf('Total Distance = %1.4f',optRoute(101)));
end
%% Calculate Function 
%calculate the mathematical distance between the points on the map using X and
%Y which is directly used as distance between cities
function d = distance(n,xy)
    comp=0;
    for i=1:length(n)-1
        if i==length(n)-1
            comp=comp+sqrt((xy(n(i),1)-xy(n(1),1))^2+(xy(n(i),2)-xy(n(1),2))^2);
        else
            comp=comp+sqrt((xy(n(i+1),1)-xy(n(i),1))^2+(xy(n(i+1),2)-xy(n(i),2))^2);
        end
    end
    d=comp;
end

%% Generate new set
%generate an individual with 100 chromosomes because there are 100 cities
function location = generatePath()
    location = randperm(100);
end

%% Cycle Cross Over
%select a random point in the array,this element will point to the same
%location in the second array, find the element in the first array that
%will point again to a vaue in the second, continue to do so until the
%start is reached, the elements of the wheel are kepts in the same
%position as the ones coresponding to the first array
%while the remaining empty locations in the array are filled with the
%elements of the second array
function C1=crossOver(parent1,parent2)

    ch1=zeros(1,length(parent1));
    wheel=[];
    ran = randi(length(parent1));
    %start with a random value in the first array parent1
    start=parent1(ran);
    wheel(end+1)=start;
    %find the position of the second element
    posSecond = find(parent1==start);
    next=parent2(posSecond);
    while next~=start
        wheel(end+1)=next;
        posi=find(parent1==next);
        next = parent2(1,posi);
    end
    %fill the remaining spaces
    for i=1:length(parent1)
           if(ismember(wheel,parent2(i))==false)
               ch1(i)=parent2(i);
           else
               add=find(wheel==parent1(i));
           
               ch1(i)=wheel(add);
           end
    end
    C1=ch1;
end

%% Order One Cross Over - work in progress - fix last row 
function ord=orderCross(parent1, parent2)

    len=length(parent1);
    firstPoint = randi([1 len]);
    secondPoint = randi([1 len]);
    if(firstPoint>secondPoint)
        mid=firstPoint;
        firstPoint=secondPoint
        secondPoint=mid
    end
    final=zeros(1,len);
    for i=1:len
        if i>=firstPoint && i<=secondPoint
            final(i)=parent1(i);
        end
    end
    keepTrackAdd=secondPoint+1;
    keepTrackArray=secondPoint+1;
    for j=1:len
        if keepTrackAdd<len
            keepTrackAdd=1;
        end
        if keepTrackArray<len
            keepTrackArray=1;
        end
        
        if ismember(parent2(keepTrackArray),final) ~= true
            final(keepTrackAdd)=parent2(keepTrackArray);
            keepTrackAdd = keepTrackAdd+1;
        end
        keepTrackArray=keepTrackArray+1;
    end
    ord=final;
end

%% Displacement Mutation
%select 2 random points, create an array with the elements inside those points, combine the elements
%outside the points, select a random gap between the curent array (the
%elements outside the points) and insert the array of the first points
%inside the current array
function M=mutate(original)
    len=length(original);
    firstPoint = randi([2 len-1]);
    secondPoint = randi([2 len-1]);
    valuesToMove=[];
    if(firstPoint>secondPoint)
        mid=firstPoint;
        firstPoint=secondPoint;
        secondPoint=mid;
    end
    for i=1:secondPoint
        if i>=firstPoint
            valuesToMove(end+1)=original(1,i);
        end
    end
    %remained elements after the cut
    afterCut=[];
    for j=1:len
        if ismember(original(1,j),valuesToMove) ~=true
            afterCut(end+1)=original(1,j);       
        end
    end
    %select gap where to insert the values
    insertHere = randi([2 length(afterCut)]);
    %concatenate the elements
    M=[afterCut(1,1:insertHere-1) valuesToMove afterCut(1,insertHere:end)];
end
%% Inversion Mutation
%take 2 random points from the total length of the individual, mirror the
%data between those points and add it back in the same location
function I=IMutate(original)
    len=length(original);
    firstPoint = randi([2 len-1]);
    secondPoint = randi([2 len-1]);
    if(firstPoint>secondPoint)
        mid=firstPoint;
        firstPoint=secondPoint;
        secondPoint=mid;
    end
    final = original;
    s=secondPoint;
    %simply swap the elements starting from the first selected point with
    %the second selected point until reach the middle
    for i=1:len
        if i>=firstPoint && i<=secondPoint
            final(1,i)=original(1,s);
            s=s-1;
        end
    end
    
    I=final;
end

%% find Best
% the function is asking for 2 parameters to see which one is the best
% measure the performance based on the distance, sort them and return the
% one with the smaller distance
function john=best(a,b,xy)
    ar = zeros(2,100);
    ar(1,1:end)=a;
    ar(2,1:end)=b;
    ar=[ar zeros(2,1)];
    for i=1:2
        ar(i,101)=distance(ar(i,:),xy);
    end
    ar=sortrows(ar,101);
    john=ar(1,1:100);
end

%% Selection of one parent
% add all the values for distance from all memebers
% generate a random number between 0 and the sum of all distances
% if the number generated is bigger than the current score and smaller than
% the score of the curent individual plus the previous individuals combined
% then select it otherwise select a random individual
% this implementation has the oposite of the roulete but selecting the
% worst becasue it will be combined with the best and teh varation will
% increase
function one = select(members,xy)
    sum=0;
    score=0;
    len=length(members);
    %one=members(randi(100),:);
    o=[];
    for i=1:len
        sum=sum+distance(members(i,:),xy);
    end
    w=round(sum);
    rando = randi(w);
    for j=1:len
        if(rando>=score && rando<=score+distance(members(j,:),xy))
            one=members(j,:);
            break;
        end
        score=distance(members(i,:),xy);
    end
    one=members(randi(100),:);
end