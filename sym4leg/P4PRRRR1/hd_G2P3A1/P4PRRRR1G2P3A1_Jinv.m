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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2020-08-07 11:02
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P4PRRRR1G2P3A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(2,1),zeros(4,3),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR1G2P3A1_Jinv: qJ has to be [3x4] (double)');
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR1G2P3A1_Jinv: xP has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4PRRRR1G2P3A1_Jinv: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR1G2P3A1_Jinv: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR1G2P3A1_Jinv: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:01:57
% EndTime: 2020-08-07 11:01:57
% DurationCPUTime: 0.15s
% Computational Cost: add. (28->20), mult. (72->60), div. (36->8), fcn. (136->26), ass. (0->56)
t28 = sin(qJ(2,4));
t19 = 0.1e1 / t28;
t29 = cos(qJ(3,4));
t59 = t19 / t29;
t35 = sin(qJ(2,3));
t21 = 0.1e1 / t35;
t40 = cos(qJ(3,3));
t58 = t21 / t40;
t37 = sin(qJ(2,2));
t22 = 0.1e1 / t37;
t41 = cos(qJ(3,2));
t57 = t22 / t41;
t39 = sin(qJ(2,1));
t23 = 0.1e1 / t39;
t42 = cos(qJ(3,1));
t56 = t23 / t42;
t55 = t28 * t29;
t54 = t35 * t40;
t53 = t37 * t41;
t52 = t39 * t42;
t51 = koppelP(1,1);
t50 = koppelP(2,1);
t49 = koppelP(3,1);
t48 = koppelP(4,1);
t47 = koppelP(1,2);
t46 = koppelP(2,2);
t45 = koppelP(3,2);
t44 = koppelP(4,2);
t43 = xP(4);
t38 = sin(qJ(3,1));
t36 = sin(qJ(3,2));
t34 = sin(qJ(3,3));
t33 = legFrame(1,2);
t32 = legFrame(2,2);
t31 = legFrame(3,2);
t30 = legFrame(4,2);
t27 = sin(qJ(3,4));
t18 = cos(t43);
t17 = sin(t43);
t16 = cos(t33);
t15 = cos(t32);
t14 = cos(t31);
t13 = cos(t30);
t12 = sin(t33);
t11 = sin(t32);
t10 = sin(t31);
t9 = sin(t30);
t8 = t12 * t38 + t16 * t52;
t7 = t11 * t36 + t15 * t53;
t6 = t10 * t34 + t14 * t54;
t5 = t12 * t52 - t16 * t38;
t4 = t11 * t53 - t15 * t36;
t3 = t10 * t54 - t14 * t34;
t2 = t13 * t55 + t9 * t27;
t1 = -t13 * t27 + t9 * t55;
t20 = [t5 * t56, t8 * t56, cos(qJ(2,1)) * t23, (-t5 * (t17 * t51 + t18 * t47) + t8 * (-t17 * t47 + t18 * t51)) * t56; t4 * t57, t7 * t57, cos(qJ(2,2)) * t22, (-t4 * (t17 * t50 + t18 * t46) + t7 * (-t17 * t46 + t18 * t50)) * t57; t3 * t58, t6 * t58, cos(qJ(2,3)) * t21, (-t3 * (t17 * t49 + t18 * t45) + t6 * (-t17 * t45 + t18 * t49)) * t58; t1 * t59, t2 * t59, cos(qJ(2,4)) * t19, (-t1 * (t17 * t48 + t18 * t44) + t2 * (-t17 * t44 + t18 * t48)) * t59;];
Jinv  = t20;
