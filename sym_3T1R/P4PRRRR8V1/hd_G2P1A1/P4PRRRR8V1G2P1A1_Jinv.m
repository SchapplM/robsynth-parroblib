% Analytische Jacobi-Matrix für parallelen Roboter
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorpose und aktiven Gelenkkoordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Input:
% xP [4x1]
%   Generalized platform coordinates
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% Jinv [4x4]
%   Analytische Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:09
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P4PRRRR8V1G2P1A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(6,1),zeros(4,3),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V1G2P1A1_Jinv: qJ has to be [3x4] (double)');
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V1G2P1A1_Jinv: xP has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P4PRRRR8V1G2P1A1_Jinv: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V1G2P1A1_Jinv: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V1G2P1A1_Jinv: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:08:53
% EndTime: 2020-08-07 11:08:53
% DurationCPUTime: 0.18s
% Computational Cost: add. (180->44), mult. (512->112), div. (20->4), fcn. (504->30), ass. (0->76)
t83 = pkin(2) * sin(qJ(3,4));
t82 = pkin(2) * cos(qJ(3,4));
t81 = pkin(2) * sin(qJ(3,3));
t80 = pkin(2) * sin(qJ(3,2));
t79 = pkin(2) * sin(qJ(3,1));
t78 = pkin(2) * cos(qJ(3,3));
t77 = pkin(2) * cos(qJ(3,2));
t76 = pkin(2) * cos(qJ(3,1));
t44 = sin(qJ(2,4));
t46 = cos(qJ(2,4));
t21 = -t46 * pkin(5) + t44 * t82;
t40 = sin(pkin(3));
t42 = cos(pkin(3));
t75 = -t21 * t42 + t40 * t83;
t52 = sin(qJ(2,3));
t58 = cos(qJ(2,3));
t23 = -t58 * pkin(5) + t52 * t78;
t74 = -t23 * t42 + t40 * t81;
t54 = sin(qJ(2,2));
t60 = cos(qJ(2,2));
t24 = -t60 * pkin(5) + t54 * t77;
t73 = -t24 * t42 + t40 * t80;
t56 = sin(qJ(2,1));
t62 = cos(qJ(2,1));
t25 = -t62 * pkin(5) + t56 * t76;
t72 = -t25 * t42 + t40 * t79;
t71 = koppelP(1,1);
t70 = koppelP(2,1);
t69 = koppelP(3,1);
t68 = koppelP(4,1);
t67 = koppelP(1,2);
t66 = koppelP(2,2);
t65 = koppelP(3,2);
t64 = koppelP(4,2);
t63 = xP(4);
t50 = legFrame(1,2);
t49 = legFrame(2,2);
t48 = legFrame(3,2);
t47 = legFrame(4,2);
t41 = cos(pkin(6));
t39 = sin(pkin(6));
t38 = cos(t63);
t37 = sin(t63);
t36 = cos(t50);
t35 = cos(t49);
t34 = cos(t48);
t33 = cos(t47);
t32 = sin(t50);
t31 = sin(t49);
t30 = sin(t48);
t29 = sin(t47);
t28 = pkin(5) * t56 + t62 * t76;
t27 = pkin(5) * t54 + t60 * t77;
t26 = pkin(5) * t52 + t58 * t78;
t22 = pkin(5) * t44 + t46 * t82;
t20 = t25 * t40 + t42 * t79;
t19 = t24 * t40 + t42 * t80;
t18 = t23 * t40 + t42 * t81;
t17 = 0.1e1 / t20;
t16 = 0.1e1 / t19;
t15 = 0.1e1 / t18;
t14 = t21 * t40 + t42 * t83;
t13 = 0.1e1 / t14;
t12 = -t39 * t28 + t72 * t41;
t11 = -t39 * t27 + t73 * t41;
t10 = -t39 * t26 + t74 * t41;
t9 = -t39 * t22 + t75 * t41;
t8 = -t12 * t36 + t32 * t20;
t7 = t12 * t32 + t36 * t20;
t6 = -t11 * t35 + t31 * t19;
t5 = t11 * t31 + t35 * t19;
t4 = -t10 * t34 + t30 * t18;
t3 = t10 * t30 + t34 * t18;
t2 = t29 * t14 - t9 * t33;
t1 = t33 * t14 + t9 * t29;
t43 = [t8 * t17, t7 * t17, (t41 * t28 + t72 * t39) * t17, (t8 * (-t37 * t71 - t38 * t67) + t7 * (-t37 * t67 + t38 * t71)) * t17; t6 * t16, t5 * t16, (t41 * t27 + t73 * t39) * t16, (t6 * (-t37 * t70 - t38 * t66) + t5 * (-t37 * t66 + t38 * t70)) * t16; t4 * t15, t3 * t15, (t41 * t26 + t74 * t39) * t15, (t4 * (-t37 * t69 - t38 * t65) + t3 * (-t37 * t65 + t38 * t69)) * t15; t2 * t13, t1 * t13, (t41 * t22 + t75 * t39) * t13, (t2 * (-t37 * t68 - t38 * t64) + t1 * (-t37 * t64 + t38 * t68)) * t13;];
Jinv  = t43;
