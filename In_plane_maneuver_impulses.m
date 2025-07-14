%INPLANE MANEUVERS
% This code is a general way of computing impulses for in-plane maneuvers
% (Hohmann, Bi-elliptic, Change of argument of perigee).
%Inputsmust be provided by user

%For Hohmann and Bi-elliptic transfer the initial orbit is the lower orbit,
%the final is the higher orbit and transfers inititate at the periapsis of
%lower orbit

% NOTE: For switch-case, use strings: "Hohmann", "Bielliptic", "Change of AOP"

%% clean
clc; close all; clear all;

%% constants
G = 6.67430e-20; % Gravitational constant (km^3/kg/s^2)
M = 5.972e24;    % Earth mass (kg)
R = 6371;     % Earth radius (km)
mu=G*M;

%% inputs
maneuver = "Change of AOP"; % "Hohmann"or "Bielliptic" or "Change of AOP" 

a_1=R+200; %semimajor axis of initial orbit
e_1=0.8; %eccentricity of first orbit
AOP_1=0; %argument of perigee of initial orbit in degrees
AOP_1 = deg2rad(AOP_1);

a_2=R+200; %semimajor axis of final orbit
e_2=0.8; %eccentricity of final orbit
AOP_2=45;%argument of perigee of final orbit
AOP_2 = deg2rad(AOP_2);

%input if bi-elliptic orbit is chosen
r_b=R+38000;% radius of B point for bi-elliptic transfer

r_p_1=a_1*(1-e_1); %radius at periapsis for initial orbit
r_a_1=a_1*(1+e_1); %radius at apoapsis for initial orbit
r_p_2=a_2*(1-e_2); %radius at periapsis for final orbit
r_a_2=a_2*(1+e_2); %radius at apoapsis for final orbit



%% validation of inputs
if a_1 <= R || a_2 <= R
    error('Semimajor axis must be greater than Earth radius');
end

if e_1 < 0 || e_1 >= 1 || e_2 < 0 || e_2 >= 1
    error('Eccentricity must be between 0 (inclusive) and 1 (exclusive)');
end

%% maneuvers 

    switch maneuver
        case "Hohmann" %Hohmann transfer and Hohmann-like transfer
            %in this specific script it was assumed that we pass from lower
            %orbit to higher orbit and that the starting point is the perigee
            %of the lower orbit

            if abs(a_1 - a_2) < 1e-6
                error('Initial and final orbits have the same semi-major axis. No transfer needed — the spacecraft is already on the target orbit.');
            end
            if AOP_1 ~=AOP_2
                error("you cannot do a purely Hohmann transfer")
            else

                a_tr=(r_p_1+r_a_2)/2; %semi major axis of transfer orbit

                V_p_tr=sqrt(mu*((2/r_p_1)-(1/a_tr))); %velocity at perigee on transfer orbit
                V_p_1=sqrt(mu*((2/r_p_1)-(1/a_1)));%velocity at perigee on initial orbit
                V_a_tr=sqrt(mu*((2/r_a_2)-(1/a_tr)));%velocity at apogee on transfer orbit
                V_a_2=sqrt(mu*((2/r_a_2)-(1/a_2)));%velocity at apogee on final orbit

                DeltaV1=V_p_tr-V_p_1;
                DeltaV2=V_a_2-V_a_tr;

                fprintf('Velocity initial orbit at perigee= %.2f km/s \n',V_p_1);
                fprintf('Velocity final orbit at apogee= %.2f km/s \n',V_a_2);
                fprintf('Transfer orbit velocity at perigee = %.2f km/s\n', V_p_tr);
                fprintf('Transfer orbit velocity at apogee = %.2f km/s\n', V_a_tr);
                fprintf('First impulse= %.2f km/s',DeltaV1);
                fprintf('Second Impulse= %.2f km/s',DeltaV2);

                %plotting the orbits in 2D
                ni = linspace(0, 2*pi, 1000); % Full orbit
                b_1=a_1*sqrt(1-(e_1)^2); %semi minor axis pf both orbits
                b_2=a_2*sqrt(1-(e_2)^2);
                e_tr=(r_a_2-r_p_1)/(r_a_2+r_p_1);%eccentricity of transfer orbit
                b_tr=a_tr*sqrt(1-(e_tr)^2); %semi minor axis of transfer orbit

                r_1=(a_1*(1-(e_1)^2))./(1+e_1*cos(ni));
                r_2=(a_2*(1-(e_2)^2))./(1+e_2*cos(ni));
                r_tr=(a_tr*(1-(e_tr)^2))./(1+e_tr*cos(ni));

                x_1=r_1.*cos(ni);%polar coordinates
                y_1=r_1.*sin(ni);

                x_2=r_2.*cos(ni);
                y_2=r_2.*sin(ni);

                x_tr=r_tr.*cos(ni);
                y_tr=r_tr.*sin(ni);

                AOP=AOP_1;

                x_rot1=x_1*cos(AOP)-y_1*sin(AOP);
                y_rot1=x_1*sin(AOP)+y_1*cos(AOP);
                x_rot2=x_2*cos(AOP)-y_2*sin(AOP);
                y_rot2=x_2*sin(AOP)+y_2*cos(AOP);
                x_rottr=x_tr*cos(AOP)-y_tr*sin(AOP);
                y_rottr=x_tr*sin(AOP)+y_tr*cos(AOP);

                figure;
                plot(x_rot1, y_rot1, 'b') % initial orbit in blue
                hold on
                plot(x_rot2, y_rot2, 'r') % final orbit in red
                plot(x_rottr,y_rottr,'g')%transfer orbit in green
                plot(0, 0, 'ko', 'MarkerFaceColor', 'k') % central body
                axis equal
                xlabel('x [km]')
                ylabel('y [km]')
                legend('Initial Orbit', 'Final Orbit','Transfer Orbit');
            end


        case "Bielliptic"
            if AOP_1 ~= AOP_2
                error("you cannot do a purely bi-elliptic transfer")
            else

                if r_b <= max(r_p_1, r_a_2)
                    error('rb must be greater than both the initial periapsis and final apoapsis radii for a valid bi-elliptic transfer.');
                else
                    a_tr_1=(r_p_1+r_b)/2;%semimajor axis first transfer orbit
                    a_tr_2=(r_b+r_p_2)/2;%semimajor axis second transfer orbit

                    V_p_tr_1=sqrt(mu*((2/r_p_1)-(1/a_tr_1))); %velocity at perigee on first transfer orbit
                    V_p_1=sqrt(mu*((2/r_p_1)-(1/a_1))); %velocity at perigee on initial orbit
                    V_a_tr_2=sqrt(mu*((2/r_b)-(1/a_tr_2))); %velocity at apogee on second transfer orbit
                    V_a_tr_1=sqrt(mu*((2/r_b)-(1/a_tr_1))); %velocity at apogee on first transfer orbit
                    V_p_2=sqrt(mu*((2/r_p_2)-(1/a_2))); %velocity at perigee on final orbit
                    V_p_tr_2=sqrt(mu*((2/r_p_2)-(1/a_tr_2))); %velocity at perigee on second transfer orbit

                    DeltaV1=V_p_tr_1-V_p_1;
                    DeltaV2=V_a_tr_2-V_a_tr_1;
                    DeltaV3=V_p_2-V_p_tr_2;

                    fprintf('Velocity initial orbit at perigee= %.2f km/s \n',V_p_1);
                    fprintf('Velocity final orbit at perigee= %.2f km/s \n',V_p_2);
                    fprintf('First impulse= %.2f km/s \n',DeltaV1);
                    fprintf('Second Impulse= %.2f km/s \n',DeltaV2);
                   if abs(DeltaV3) < 1e-6
                       fprintf('Third Impulse= 0.00 km/s (no burn needed, orbits already match)\n');
                   else
                       fprintf('Third Impulse= %.2f km/s \n', DeltaV3);
                   end

                    %plotting the orbits in 2D
                    ni = linspace(0, 2*pi, 1000); % Full orbit
                    b_1=a_1*sqrt(1-(e_1)^2); %semi minor axis of both orbits
                    b_2=a_2*sqrt(1-(e_2)^2);
                    e_tr_1=(r_b-r_p_1)/(r_b+r_p_1);%eccentricity of first transfer orbit
                    e_tr_2=(r_b-r_p_2)/(r_b+r_p_2);%eccentricity of second transfer orbit
                    b_tr_1=a_tr_1*sqrt(1-(e_tr_1)^2); %semi minor axis of first transfer orbit
                    b_tr_2=a_tr_2*sqrt(1-(e_tr_2)^2); %semi minor axis of second transfer orbit

                    r_1=(a_1*(1-(e_1)^2))./(1+e_1*cos(ni));
                    r_2=(a_2*(1-(e_2)^2))./(1+e_2*cos(ni));
                    r_tr_1=(a_tr_1*(1-(e_tr_1)^2))./(1+e_tr_1*cos(ni));
                    r_tr_2=(a_tr_2*(1-(e_tr_2)^2))./(1+e_tr_2*cos(ni));

                    x_1=r_1.*cos(ni);%polar coordinates
                    y_1=r_1.*sin(ni);

                    x_2=r_2.*cos(ni);
                    y_2=r_2.*sin(ni);

                    x_tr_1=r_tr_1.*cos(ni);
                    y_tr_1=r_tr_1.*sin(ni);
                    x_tr_2=r_tr_2.*cos(ni);
                    y_tr_2=r_tr_2.*sin(ni);

                    AOP=AOP_1;

                    x_rot1=x_1*cos(AOP)-y_1*sin(AOP);
                    y_rot1=x_1*sin(AOP)+y_1*cos(AOP);
                    x_rot2=x_2*cos(AOP)-y_2*sin(AOP);
                    y_rot2=x_2*sin(AOP)+y_2*cos(AOP);
                    x_rottr1=x_tr_1*cos(AOP)-y_tr_1*sin(AOP);
                    y_rottr1=x_tr_1*sin(AOP)+y_tr_1*cos(AOP);
                    x_rottr2=x_tr_2*cos(AOP)-y_tr_2*sin(AOP);
                    y_rottr2=x_tr_2*sin(AOP)+y_tr_2*cos(AOP);


                    figure;
                    plot(x_rot1, y_rot1, 'b') % initial orbit in blue
                    hold on
                    plot(x_rot2, y_rot2, 'r') % final orbit in red
                    plot(x_rottr1,y_rottr1,'g')%first transfer orbit in green
                    plot(x_rottr2,y_rottr2,'m')%second transfer orbit in magenta
                    plot(0, 0, 'ko', 'MarkerFaceColor', 'k') % central body
                    axis equal
                    xlabel('x [km]')
                    ylabel('y [km]')
                    legend('Initial Orbit', 'Final Orbit','First Transfer Orbit','Second Transfer Orbit');

                    if e_1==0 && e_2==0
                        if r_p_2/r_p_1<11.94
                            disp("Hohmann transfer is more efficient")
                        end
                    end
                end


            end

        case "Change of AOP"
            if (e_2~=e_1) || (a_1~=a_2)
                error("this is not a change of AOP, eccentricity and semimajor axis should be the same for both initial and final orbit")
            else
                DAOP=AOP_2-AOP_1;
                e=e_1;
                a=a_2;
                p=a*(1-e^2);
                DeltaV=2*e*(sqrt(mu/p))*sin(DAOP/2);
                fprintf('Impulse= %.2f km/s',DeltaV);

                %plotting the orbits in 2D
                ni = linspace(0, 2*pi, 1000); % Full orbit
                b=a*sqrt(1-(e)^2); %semi minor axis pf both orbits
                r=(a*(1-(e)^2))./(1+e_1*cos(ni));
                rotate_orbit = @(r, ni, AOP) [r .* cos(ni + AOP); r .* sin(ni + AOP)];

                % Initial orbit
                r1 = rotate_orbit(r, ni, AOP_1);

                % Final orbit
                r2 = rotate_orbit(r, ni,AOP_2);

                % Plot
                figure; 
                hold on; axis equal; grid on;
                plot(r1(1,:), r1(2,:), 'b', 'LineWidth', 2);
                plot(r2(1,:), r2(2,:), 'r', 'LineWidth', 2);
                plot(0,0,'ko','MarkerFaceColor','k');
                

                legend('Initial Orbit', 'Final Orbit', 'Location', 'best');
                title('Change of Argument of Perigee Maneuver');
                xlabel('x [km]'); ylabel('y [km]');
            end



    end











    



