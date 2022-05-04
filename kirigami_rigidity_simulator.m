% Interactive kirigami rigidity simulator
%
% Usage:
% left click: rigidify a tile/hole
% right click: soften a tile/hole 
% - does not work if the tile/hole rigidity is already determined by some
%   other tiles and holes

global pt clicktype

% Pattern size
m = 6;
n = 6;

[vertices,holes] = generate_kirigami(m,n);

scale = 100;
A = zeros(m,n); % tile rigidity
B = zeros(m-1,n-1); % hole rigidity

fprintf('Usage:\n');
fprintf('- Left click: rigidify a tile/hole\n');
fprintf('- Right click: soften a tile/hole\n');

last_chosen = [];
f=figure(1);
while (sum(sum(A)) ~= m*n) || (sum(sum(B)) ~= (m-1)*(n-1))
    %% click a tile/hole in the figure to rigidify it
    clf(1);
    set(f,'WindowButtonDownFcn',@detectpoint);
%     imshow(im);
    plot_kirigami(vertices*scale,holes*scale,A,B,last_chosen);
    waitforbuttonpress;
    
    % Rigidify the corresponding quad or hole based on the clicked coordinates
    k = 1;
    check = 0;
    check2 = 0;
    while check == 0 && k <= length(vertices)
         check = inpolygon(pt(1,1),pt(1,2),vertices(k,1:2:end)*scale,vertices(k,2:2:end)*scale);
        if check == 0
            k = k+1;
        else
            check2 = 1;
        end
    end
    
    k2 = 1;
    while check2 == 0 && k2 <= length(holes)
         check2 = inpolygon(pt(1,1),pt(1,2),holes(k2,1:2:end)*scale,holes(k2,2:2:end)*scale);
        if check2 == 0
            k2 = k2+1;
        else
            check2 = 1;
        end
    end
        
    if check == 1
        i = floor((k-1)/size(A,2))+1;
        j = k - (i-1)*size(A,2);
        if strcmpi(clicktype,'normal')
            A(i,j) = 1;
        else
            A(i,j) = 0;
        end
        last_chosen = [1,k];
        
    elseif check2 == 1
        i = floor((k2-1)/size(B,2))+1;
        j = k2 - (i-1)*size(B,2);
        if strcmpi(clicktype,'normal')
            B(i,j) = 1;
        else
            B(i,j) = 0;
        end
        last_chosen = [2,k2];
    end
    
    
    % Check and update the overall tile and hole rigidity
    % A tile/hole becomes rigid if one of its angles is uniquely determined
    A_old = A;
    B_old = B;
    converge = 0;
    while converge == 0
        
        A_temp = zeros(m+2,n+2);
        B_temp = zeros(m+1,n+1);
        A_temp(2:end-1,2:end-1) = A;
        B_temp(2:end-1,2:end-1) = B;
        
        % check if a hole becomes rigid
        for i = 1:m-1
            for j = 1:n-1
%                 B(i,j) = max([B(i,j),...
%                     B(i-1,j)*A(i,j)*A(i,j+1),...
%                     B(i+1,j)*A(i+1,j)*A(i+1,j+1),...
%                     B(i,j-1)*A(i,j)*A(i+1,j),...
%                     B(i,j+1)*A(i,j+1)*A(i+1,j+1)]);
                B(i,j) = max([B_temp(i+1,j+1),...
                    B_temp(i,j+1)*A_temp(i+1,j+1)*A_temp(i+1,j+1+1),...
                    B_temp(i+1+1,j+1)*A_temp(i+1+1,j+1)*A_temp(i+1+1,j+1+1),...
                    B_temp(i+1,j-1+1)*A_temp(i+1,j+1)*A_temp(i+1+1,j+1),...
                    B_temp(i+1,j+1+1)*A_temp(i+1,j+1+1)*A_temp(i+1+1,j+1+1)]);
            end
        end
        
        % check if a quad becomes rigid
        for i = 1:m
            for j = 1:n
%                 A(i,j) = max([A(i,j),...
%                     A(i-1,j)*B(i-1,j-1)*B(i-1,j),...
%                     A(i,j+1)*B(i,j)*B(i-1,j),...
%                     A(i+1,j)*B(i,j)*B(i,j-1),...
%                     A(i,j-1)*B(i-1,j-1)*B(i,j-1)]);
                A(i,j) = max([A_temp(i+1,j+1),...
                    A_temp(i-1+1,j+1)*B_temp(i-1+1,j-1+1)*B_temp(i-1+1,j+1),...
                    A_temp(i+1,j+1+1)*B_temp(i+1,j+1)*B_temp(i-1+1,j+1),...
                    A_temp(i+1+1,j+1)*B_temp(i+1,j+1)*B_temp(i+1,j-1+1),...
                    A_temp(i+1,j-1+1)*B_temp(i-1+1,j-1+1)*B_temp(i+1,j-1+1)]);
            end
        end
        
        if norm(A-A_old,2) + norm(B-B_old) == 0
            converge = 1;
        else 
           A_old = A;
           B_old = B;
        end
    end
end

clf(1);
plot_kirigami(vertices*scale,holes*scale,A,B,last_chosen);

fprintf('Completely rigidified!\n');

%%
function detectpoint(h,~)
    global pt clicktype
    clicktype = get(gcf,'SelectionType');
    pt = get(gca,'CurrentPoint');
    fprintf([clicktype, ' click: %d %d\n'], pt(1,1),pt(1,2));
end


function [vertices,holes] = generate_kirigami(m,n)
    vertices = [];
    theta = pi/12;
    temp = 1;
    for i = 1:m
        for j = 1:n
            % square id: the (i,j)-th square = (i-1)*m+j
            % 4 - 3
            % |   |
            % 1 - 2

            temp_quads = [j-1,i-1, j,i-1, j,i, j-1,i];

            % shift the quads to center at origin
            temp_quads(:,1:2:end) = temp_quads(:,1:2:end) - (j-1/2);
            temp_quads(:,2:2:end) = temp_quads(:,2:2:end) - (i-1/2);
            % rotate the quads by theta or -theta
            if mod(i+j,2) == 1
                temp_quads_rotated1 = cos(theta)*temp_quads(:,1:2:end)-sin(theta)*temp_quads(:,2:2:end);
                temp_quads_rotated2 = sin(theta)*temp_quads(:,1:2:end)+cos(theta)*temp_quads(:,2:2:end);
            else
                temp_quads_rotated1 = cos(-theta)*temp_quads(:,1:2:end)-sin(-theta)*temp_quads(:,2:2:end);
                temp_quads_rotated2 = sin(-theta)*temp_quads(:,1:2:end)+cos(-theta)*temp_quads(:,2:2:end);

            end
            temp_quads_rotated = zeros(size(temp_quads));
            temp_quads_rotated(:,1:2:end) = temp_quads_rotated1;
            temp_quads_rotated(:,2:2:end) = temp_quads_rotated2;

            % rescale the quads so that they are contained in the unit square
            temp_quads_rotated = temp_quads_rotated/(sin(theta)+cos(theta));

            % shift them back to the desired position
            temp_quads_rotated(:,1:2:end) = temp_quads_rotated(:,1:2:end) + (j-1/2);
            temp_quads_rotated(:,2:2:end) = temp_quads_rotated(:,2:2:end) + (i-1/2);

            vertices(temp,:) = temp_quads_rotated;
            temp = temp+1;
        end
    end
    
    holes = [];
    for i = 1:m-1
        for j = 1:n-1
            k = (i-1)*n+j;
            k2 = k + n+1;
            v = [vertices(k,1:2:end)', vertices(k,2:2:end)'];
            v2 = [vertices(k2,1:2:end)', vertices(k2,2:2:end)'];
            
            if mod(i+j,2) == 0
                holes = [holes; v(3,1:2),v2(2,1:2),v2(1,1:2),v(4,1:2)];
            else
                holes = [holes; v(2,1:2),v2(1,1:2),v2(4,1:2),v(3,1:2)];
            end
        end
    end
    
end

function plot_kirigami(vertices,holes,A,B,last_chosen)

    axis equal on;
    hold on;
    for k = 1:length(vertices)
        v = [vertices(k,1:2:end)', vertices(k,2:2:end)'];
%         (i-1)*m+j
        i = floor((k-1)/size(A,2))+1;
        j = k - (i-1)*size(A,2);
        if A(i,j)
            fill(v(:,1), v(:,2), [255 230 213]/255, 'LineWidth',1);
        else
            fill(v(:,1), v(:,2), [255 255 255]/255, 'LineWidth',1);
        end
        if i ~= size(A,1) && j ~= size(A,2)
            k2 = (i-1)*size(B,2)+j;
            if B(i,j)
                fill(holes(k2,1:2:end),holes(k2,2:2:end),[180 252 171]/255, 'LineWidth',1);
            end
        end
    end
    if ~isempty(last_chosen)
        k = last_chosen(2);
        switch last_chosen(1)
            case 1
                % quad
                plot(vertices(k,[1:2:end,1]),vertices(k,[2:2:end,2]),'-','Color',[201 0 22]/255,'LineWidth',2);
            case 2
                % negative space
                plot(holes(k,[1:2:end,1]),holes(k,[2:2:end,2]),'-','Color',[201 0 22]/255,'LineWidth',2);
        end
    end
end