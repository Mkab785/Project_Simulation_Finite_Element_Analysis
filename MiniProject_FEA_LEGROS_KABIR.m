clear;
clc;

%% 1. Data Input
Nr = input('Enter the number of bars (Nr): ');
Nn = input('Enter the number of nodes (Nn): ');

Module = zeros(Nr, 1);
Section = zeros(Nr, 1);
Coord = zeros(Nn, 2);
Connec = zeros(Nr, 2);

disp(' ');
disp('--- Entering Bar Properties (E, A) ---');
for i = 1:Nr
    Module(i) = input(sprintf('Enter Young''s Modulus (E) for bar %d: ', i));
    Section(i) = input(sprintf('Enter Cross-sectional Area (A) for bar %d: ', i));
end

disp(' ');
disp('--- Entering Coordinates ---');
coord_unique = {};
for i = 1:Nn
    unique_coord = false;
    while ~unique_coord
        x = input(sprintf('Enter the ''x'' coordinate for node %d: ', i));
        y = input(sprintf('Enter the ''y'' coordinate for node %d: ', i));
        
        point_actuel = [x, y];
        is_duplicate = false;
        
        for k = 1:length(coord_unique)
            if isequal(coord_unique{k}, point_actuel)
                is_duplicate = true;
                break;
            end
        end
        
        if is_duplicate
            disp(sprintf('Error: Point (%g, %g) already exists. Enter different coordinates.', x, y));
        else
            unique_coord = true;q
            Coord(i, 1) = x;
            Coord(i, 2) = y;
            coord_unique{end+1} = point_actuel;
        end
    end
end

disp(' ');
disp('--- Entering Connectivity ---');
barre_unique = {};
for i = 1:Nr
    pt_depart = 0;
    pt_arrivee = 0;
    unique_bar = false;
    
    while ~unique_bar || pt_depart == 0 || pt_depart == pt_arrivee
        
        while ~(pt_depart > 0 && pt_depart <= Nn)
            pt_depart = input(sprintf('Enter the start node number (1 to %d) for bar %d: ', Nn, i));
        end
        
        while ~(pt_arrivee > 0 && pt_arrivee <= Nn)
            pt_arrivee = input(sprintf('Enter the end node number (1 to %d) for bar %d: ', Nn, i));
        end
        
        if pt_depart == pt_arrivee
            disp('Error: Start and end nodes must be different.');
            pt_depart = 0;
            pt_arrivee = 0;
            continue;
        end
        
        barre_actuelle = sort([pt_depart, pt_arrivee]);
        is_duplicate = false;
        
        for k = 1:length(barre_unique)
            if isequal(barre_unique{k}, barre_actuelle)
                is_duplicate = true;
                break;
            end
        end
        
        if is_duplicate
            disp('This bar already exists.');
            pt_depart = 0;
            pt_arrivee = 0;
            continue;
        else
            unique_bar = true;
        end
    end
    
    barre_unique{end+1} = barre_actuelle;
    Connec(i, 1) = pt_depart;
    Connec(i, 2) = pt_arrivee;
end

fprintf('\nNumber of bars (Nr): %d\n', Nr);
fprintf('Number of nodes (Nn): %d\n', Nn);
disp('Young''s Modulus (Module):'); disp(Module);
disp('Cross-sectional Area (Section):'); disp(Section);
disp('Nodal Coordinates (Coord):'); disp(Coord);
disp('Connectivity (Connec):'); disp(Connec);

%% 2. Length and Orientation
Length = zeros(Nr, 1);
Orient = zeros(Nr, 1);

for i = 1:Nr
    node_start = Connec(i, 1);
    node_end = Connec(i, 2);
    
    x1 = Coord(node_start, 1);
    y1 = Coord(node_start, 2);
    x2 = Coord(node_end, 1);
    y2 = Coord(node_end, 2);
    
    L = sqrt((x2 - x1)^2 + (y2 - y1)^2);
    Length(i) = L;
    
    theta = atan2(y2 - y1, x2 - x1);
    Orient(i) = theta;
end

disp(' ');
disp('Lengths:'); disp(Length');
disp('Orientations (in radians):'); disp(Orient');
disp('Orientations (in degrees):'); disp(rad2deg(Orient)');

%% 3. Global Stiffness Matrix
Dof = 2 * Nn;
K = zeros(Dof, Dof);

for i = 1:Nr
    E = Module(i);
    A = Section(i);
    L = Length(i);
    theta = Orient(i);
    
    c = cos(theta);
    s = sin(theta);
    
    Ke_prime = (E * A / L) * [
         c*c,  c*s, -c*c, -c*s;
         c*s,  s*s, -c*s, -s*s;
        -c*c, -c*s,  c*c,  c*s;
        -c*s, -s*s,  c*s,  s*s
    ];

    fprintf('\nStiffness Matrix for bar %d:\n', i);
    disp(Ke_prime);
    
    n1 = Connec(i, 1);
    n2 = Connec(i, 2);
    
    dof_map = [2*n1 - 1, 2*n1, 2*n2 - 1, 2*n2];
    
    K(dof_map, dof_map) = K(dof_map, dof_map) + Ke_prime;
end

disp(' ');
disp('Global Stiffness Matrix (K):');
disp(K);

%% 4. Boundary Conditions and Force Vector 
BC_Ux = zeros(Nn, 1);
BC_Uy = zeros(Nn, 1);
VAL_Ux = zeros(Nn, 1);
VAL_Uy = zeros(Nn, 1);
VAL_Fx = zeros(Nn, 1);
VAL_Fy = zeros(Nn, 1);
F = zeros(Dof, 1);

disp(' ');
disp('--- Entering Boundary Conditions and Loads ---');

for i = 1:Nn
    fprintf('\n--- Node %d ---\n', i);

    while true
        choix_ux = input(sprintf('Is U_x%d constrained/blocked? (0 = No and 1 = Yes): ', i));
        if choix_ux == 0 || choix_ux == 1
            break;
        else
            disp('Invalid input. Enter 0 or 1.');
        end
    end
    
    if choix_ux == 1
        BC_Ux(i) = 1;
        VAL_Ux(i) = input(sprintf('Enter the value of U_x%d (e.g., 0): ', i));
    else
        BC_Ux(i) = 0;
        
        while true
            choix_fx = input(sprintf('Is a Force F_x%d applied? (0 = No and 1 = Yes): ', i));
            if choix_fx == 0 || choix_fx == 1
                break;
            else
                disp('Invalid input. Enter 0 or 1.');
            end
        end
        if choix_fx == 1
            VAL_Fx(i) = input(sprintf('Enter the value of F_x%d (in Newtons): ', i));
            F(2*i - 1) = VAL_Fx(i);
        end
    end

    while true
        choix_uy = input(sprintf('Is U_y%d constrained/blocked? (0 = No and 1 = Yes): ', i));
        if choix_uy == 0 || choix_uy == 1
            break;
        else
            disp('Invalid input. Enter 0 or 1.');
        end
    end
    
    if choix_uy == 1
        BC_Uy(i) = 1;
        VAL_Uy(i) = input(sprintf('Enter the value of U_y%d (e.g., 0): ', i));
    else
        BC_Uy(i) = 0;
        
        while true
            choix_fy = input(sprintf('Is a Force F_y%d applied? (0 = No and 1 = Yes): ', i));
            if choix_fy == 0 || choix_fy == 1
                break;
            else
                disp('Invalid input. Enter 0 or 1.');
            end
        end
        if choix_fy == 1
            VAL_Fy(i) = input(sprintf('Enter the value of F_y%d (in Newtons): ', i));
            F(2*i) = VAL_Fy(i);
        end
    end
end

Row_Column_elimin = zeros(Dof, 1);
N_BC_U = 0;
for i = 1:Nn
    if BC_Ux(i) == 1
        Row_Column_elimin(2*i - 1) = 1;
        N_BC_U = N_BC_U + 1;
    end
    if BC_Uy(i) == 1
        Row_Column_elimin(2*i) = 1;
        N_BC_U = N_BC_U + 1;
    end
end

disp(' ');
disp('Elimination Map (1=Fixed DOF):'); disp(Row_Column_elimin');
fprintf('Number of Fixed DOFs (N_BC_U): %d\n', N_BC_U);
disp('Full Force Vector (F):'); disp(F');

free_dofs = find(Row_Column_elimin == 0);

K_red = K(free_dofs, free_dofs);
F_red = F(free_dofs);

disp(' ');
disp('Reduced Stiffness Matrix (K_red):'); disp(K_red);
disp('Reduced Force Vector (F_red):'); disp(F_red');

%% 5. Solve for Displacements
q = zeros(Dof, 1);
q_red = K_red \ F_red;

k = 1;
for i = 1:Dof
    if Row_Column_elimin(i) == 0
        q(i) = q_red(k);
        k = k + 1;
    else
        if mod(i, 2) ~= 0
            q(i) = VAL_Ux((i+1)/2);
        else
            q(i) = VAL_Uy(i/2);
        end
    end
end

disp(' ');
disp('Reduced Nodal Displacement Vector (q_red):'); disp(q_red');
disp('Final Nodal Displacement Vector (q):'); disp(q');

%% 6. Visualization
echelle = 10000; 

figure;
hold on;
title(sprintf('Deformed Structure (Scale x%d) - | Initial Structure --', echelle));

for i = 1:Nr
    n1 = Connec(i, 1);
    n2 = Connec(i, 2);
    
    plot([Coord(n1, 1), Coord(n2, 1)], [Coord(n1, 2), Coord(n2, 2)], 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    
    x1 = Coord(n1, 1) + q(2*n1 - 1) * echelle;
    y1 = Coord(n1, 2) + q(2*n1) * echelle;
    x2 = Coord(n2, 1) + q(2*n2 - 1) * echelle;
    y2 = Coord(n2, 2) + q(2*n2) * echelle;
    
    plot([x1, x2], [y1, y2], 'r-', 'LineWidth', 1.5);
end

for i = 1:Nn
    x = Coord(i, 1) + q(2*i - 1) * echelle;
    y = Coord(i, 2) + q(2*i) * echelle;
    
    plot(x, y, 'ro', 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
    text(x + 0.1, y + 0.1, num2str(i), 'Color', 'red', 'FontSize', 12, 'FontWeight', 'bold', 'HandleVisibility', 'off');
end

axis equal;
grid on;
xlabel('X Coordinate');
ylabel('Y Coordinate');
hold off;