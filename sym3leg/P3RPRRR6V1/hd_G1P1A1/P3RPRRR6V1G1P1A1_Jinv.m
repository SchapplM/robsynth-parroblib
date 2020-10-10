% Analytische Jacobi-Matrix für parallelen Roboter
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorpose und aktiven Gelenkkoordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% Jinv [3x3]
%   Analytische Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:35
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RPRRR6V1G1P1A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(7,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G1P1A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G1P1A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G1P1A1_Jinv: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G1P1A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G1P1A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:35:15
% EndTime: 2020-08-06 18:35:15
% DurationCPUTime: 0.25s
% Computational Cost: add. (234->82), mult. (174->55), div. (6->3), fcn. (78->54), ass. (0->42)
t62 = 2 * pkin(2);
t41 = (pkin(6) + pkin(5));
t61 = -2 * t41;
t59 = 2 * t41;
t47 = (qJ(1,3) + legFrame(3,3));
t28 = (pkin(7) + t47);
t17 = qJ(3,3) + t28;
t18 = -qJ(3,3) + t28;
t58 = sin(t17) + sin(t18);
t48 = (qJ(1,2) + legFrame(2,3));
t29 = (pkin(7) + t48);
t21 = qJ(3,2) + t29;
t22 = -qJ(3,2) + t29;
t57 = sin(t21) + sin(t22);
t49 = (qJ(1,1) + legFrame(1,3));
t30 = (pkin(7) + t49);
t25 = qJ(3,1) + t30;
t26 = -qJ(3,1) + t30;
t56 = sin(t25) + sin(t26);
t42 = 2 * qJ(3,3);
t55 = 0.1e1 / (pkin(3) * sin(t42) + sin(qJ(3,3)) * t62 + (sin((pkin(7) + qJ(3,3))) + sin((-pkin(7) + qJ(3,3)))) * pkin(1)) / 0.2e1;
t43 = 2 * qJ(3,2);
t54 = 0.1e1 / (pkin(3) * sin(t43) + sin(qJ(3,2)) * t62 + (sin((pkin(7) + qJ(3,2))) + sin((-pkin(7) + qJ(3,2)))) * pkin(1)) / 0.2e1;
t44 = 2 * qJ(3,1);
t53 = 0.1e1 / (pkin(3) * sin(t44) + sin(qJ(3,1)) * t62 + (sin((pkin(7) + qJ(3,1))) + sin((-pkin(7) + qJ(3,1)))) * pkin(1)) / 0.2e1;
t52 = cos(t17) + cos(t18);
t51 = cos(t21) + cos(t22);
t50 = cos(t25) + cos(t26);
t46 = 0.2e1 * pkin(1);
t36 = -qJ(3,1) + t49;
t35 = qJ(3,1) + t49;
t34 = -qJ(3,2) + t48;
t33 = qJ(3,2) + t48;
t32 = -qJ(3,3) + t47;
t31 = qJ(3,3) + t47;
t27 = -2 * qJ(3,1) + t30;
t24 = t44 + t30;
t23 = -2 * qJ(3,2) + t29;
t20 = t43 + t29;
t19 = -2 * qJ(3,3) + t28;
t16 = t42 + t28;
t1 = [(t50 * t62 + (cos(t36) + cos(t35)) * t46 + t56 * t59 + (0.2e1 * cos(t30) + cos(t27) + cos(t24)) * pkin(3)) * t53, (t56 * t62 + (sin(t36) + sin(t35)) * t46 + t50 * t61 + (sin(t27) + sin(t24) + 0.2e1 * sin(t30)) * pkin(3)) * t53, 1; (t51 * t62 + (cos(t34) + cos(t33)) * t46 + t57 * t59 + (cos(t23) + cos(t20) + 0.2e1 * cos(t29)) * pkin(3)) * t54, (t57 * t62 + (sin(t34) + sin(t33)) * t46 + t51 * t61 + (sin(t23) + sin(t20) + 0.2e1 * sin(t29)) * pkin(3)) * t54, 1; (t52 * t62 + (cos(t32) + cos(t31)) * t46 + t58 * t59 + (cos(t19) + cos(t16) + 0.2e1 * cos(t28)) * pkin(3)) * t55, (t58 * t62 + (sin(t32) + sin(t31)) * t46 + t52 * t61 + (sin(t19) + sin(t16) + 0.2e1 * sin(t28)) * pkin(3)) * t55, 1;];
Jinv  = t1;
