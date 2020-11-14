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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
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
% Datum: 2020-08-06 17:30
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3PRRRR8V1G4P2A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(6,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G4P2A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G4P2A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G4P2A1_Jinv: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G4P2A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G4P2A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:29:55
% EndTime: 2020-08-06 17:29:55
% DurationCPUTime: 0.15s
% Computational Cost: add. (180->45), mult. (510->105), div. (9->3), fcn. (519->34), ass. (0->71)
t79 = pkin(2) * sin(qJ(3,3));
t78 = pkin(2) * sin(qJ(3,2));
t77 = pkin(2) * sin(qJ(3,1));
t76 = pkin(2) * cos(qJ(3,3));
t75 = pkin(2) * cos(qJ(3,2));
t74 = pkin(2) * cos(qJ(3,1));
t60 = sin(qJ(2,3));
t66 = cos(qJ(2,3));
t22 = -t66 * pkin(5) + t60 * t76;
t47 = sin(pkin(3));
t49 = cos(pkin(3));
t73 = -t22 * t49 + t47 * t79;
t62 = sin(qJ(2,2));
t68 = cos(qJ(2,2));
t23 = -t68 * pkin(5) + t62 * t75;
t72 = -t23 * t49 + t47 * t78;
t64 = sin(qJ(2,1));
t70 = cos(qJ(2,1));
t24 = -t70 * pkin(5) + t64 * t74;
t71 = -t24 * t49 + t47 * t77;
t58 = legFrame(1,2);
t57 = legFrame(2,2);
t56 = legFrame(3,2);
t55 = legFrame(1,1);
t54 = legFrame(2,1);
t53 = legFrame(3,1);
t52 = legFrame(1,3);
t51 = legFrame(2,3);
t50 = legFrame(3,3);
t48 = cos(pkin(6));
t46 = sin(pkin(6));
t45 = cos(t58);
t44 = cos(t57);
t43 = cos(t56);
t42 = sin(t58);
t41 = sin(t57);
t40 = sin(t56);
t39 = cos(t55);
t38 = cos(t54);
t37 = cos(t53);
t36 = cos(t52);
t35 = cos(t51);
t34 = cos(t50);
t33 = sin(t55);
t32 = sin(t54);
t31 = sin(t53);
t30 = sin(t52);
t29 = sin(t51);
t28 = sin(t50);
t27 = pkin(5) * t64 + t70 * t74;
t26 = pkin(5) * t62 + t68 * t75;
t25 = pkin(5) * t60 + t66 * t76;
t18 = t24 * t47 + t49 * t77;
t17 = t23 * t47 + t49 * t78;
t16 = t22 * t47 + t49 * t79;
t15 = 0.1e1 / t18;
t14 = 0.1e1 / t17;
t13 = 0.1e1 / t16;
t12 = -t46 * t27 + t71 * t48;
t11 = -t46 * t26 + t72 * t48;
t10 = -t46 * t25 + t73 * t48;
t9 = t27 * t48 + t71 * t46;
t8 = t26 * t48 + t72 * t46;
t7 = t25 * t48 + t73 * t46;
t6 = t12 * t36 - t9 * t30;
t5 = t11 * t35 - t8 * t29;
t4 = t10 * t34 - t7 * t28;
t3 = -t45 * t18 + (t12 * t30 + t36 * t9) * t42;
t2 = -t44 * t17 + (t11 * t29 + t35 * t8) * t41;
t1 = -t43 * t16 + (t10 * t28 + t34 * t7) * t40;
t19 = [(((-t30 * t46 + t48 * t36) * t27 + t71 * (t30 * t48 + t36 * t46)) * t45 + t42 * t18) * t15, (t3 * t33 - t39 * t6) * t15, (-t3 * t39 - t33 * t6) * t15; (((-t29 * t46 + t48 * t35) * t26 + t72 * (t29 * t48 + t35 * t46)) * t44 + t41 * t17) * t14, (t2 * t32 - t38 * t5) * t14, (-t2 * t38 - t32 * t5) * t14; (((-t28 * t46 + t48 * t34) * t25 + t73 * (t28 * t48 + t34 * t46)) * t43 + t40 * t16) * t13, (t1 * t31 - t37 * t4) * t13, (-t1 * t37 - t31 * t4) * t13;];
Jinv  = t19;
