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
%   pkin=[a2,a4,alpha2,d1,d2,d4]';
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
% Datum: 2020-08-06 22:53
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRRRR10V1G2P2A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(6,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR10V1G2P2A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR10V1G2P2A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRRRR10V1G2P2A1_Jinv: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR10V1G2P2A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR10V1G2P2A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 22:52:45
% EndTime: 2020-08-06 22:52:45
% DurationCPUTime: 0.16s
% Computational Cost: add. (90->42), mult. (279->102), div. (9->3), fcn. (243->26), ass. (0->60)
t24 = sin(qJ(3,3));
t59 = pkin(2) * t24;
t27 = sin(qJ(3,2));
t58 = pkin(2) * t27;
t30 = sin(qJ(3,1));
t57 = pkin(2) * t30;
t25 = sin(qJ(2,3));
t56 = t25 * pkin(6);
t28 = sin(qJ(2,2));
t55 = t28 * pkin(6);
t31 = sin(qJ(2,1));
t54 = t31 * pkin(6);
t20 = cos(pkin(3));
t33 = cos(qJ(3,3));
t34 = cos(qJ(2,3));
t7 = -t25 * t33 * pkin(2) + t34 * pkin(6);
t53 = t7 * t20;
t36 = cos(qJ(3,2));
t37 = cos(qJ(2,2));
t8 = -t28 * t36 * pkin(2) + t37 * pkin(6);
t52 = t8 * t20;
t39 = cos(qJ(3,1));
t40 = cos(qJ(2,1));
t9 = -t31 * t39 * pkin(2) + t40 * pkin(6);
t51 = t9 * t20;
t50 = t33 * t34;
t49 = t36 * t37;
t48 = t39 * t40;
t35 = cos(qJ(1,3));
t47 = t35 * t59;
t38 = cos(qJ(1,2));
t46 = t38 * t58;
t41 = cos(qJ(1,1));
t45 = t41 * t57;
t44 = t20 * t59;
t43 = t20 * t58;
t42 = t20 * t57;
t32 = sin(qJ(1,1));
t29 = sin(qJ(1,2));
t26 = sin(qJ(1,3));
t23 = legFrame(1,2);
t22 = legFrame(2,2);
t21 = legFrame(3,2);
t19 = sin(pkin(3));
t18 = cos(t23);
t17 = cos(t22);
t16 = cos(t21);
t15 = sin(t23);
t14 = sin(t22);
t13 = sin(t21);
t12 = pkin(2) * t48 + t54;
t11 = pkin(2) * t49 + t55;
t10 = pkin(2) * t50 + t56;
t6 = -t32 * t12 + t41 * t51;
t5 = -t29 * t11 + t38 * t52;
t4 = -t26 * t10 + t35 * t53;
t3 = 0.1e1 / (pkin(1) * t51 + (-pkin(5) * t54 + (pkin(1) * t30 - pkin(5) * t48) * pkin(2)) * t19);
t2 = 0.1e1 / (pkin(1) * t52 + (-pkin(5) * t55 + (pkin(1) * t27 - pkin(5) * t49) * pkin(2)) * t19);
t1 = 0.1e1 / (pkin(1) * t53 + (-pkin(5) * t56 + (pkin(1) * t24 - pkin(5) * t50) * pkin(2)) * t19);
t60 = [((t15 * t9 + t18 * t45) * t19 + t6 * t18 - t15 * t42) * t3, ((-t15 * t45 + t18 * t9) * t19 - t6 * t15 - t18 * t42) * t3, (-t41 * t12 + (-t19 * t57 - t51) * t32) * t3; ((t14 * t8 + t17 * t46) * t19 + t5 * t17 - t14 * t43) * t2, ((-t14 * t46 + t17 * t8) * t19 - t5 * t14 - t17 * t43) * t2, (-t38 * t11 + (-t19 * t58 - t52) * t29) * t2; ((t13 * t7 + t16 * t47) * t19 + t4 * t16 - t13 * t44) * t1, ((-t13 * t47 + t16 * t7) * t19 - t4 * t13 - t16 * t44) * t1, (-t35 * t10 + (-t19 * t59 - t53) * t26) * t1;];
Jinv  = t60;
