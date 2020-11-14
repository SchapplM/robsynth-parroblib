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
% Datum: 2020-08-06 22:13
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRRRR10V1G1P1A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(6,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR10V1G1P1A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR10V1G1P1A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRRRR10V1G1P1A1_Jinv: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR10V1G1P1A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR10V1G1P1A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 22:12:24
% EndTime: 2020-08-06 22:12:24
% DurationCPUTime: 0.15s
% Computational Cost: add. (90->45), mult. (261->97), div. (9->3), fcn. (249->26), ass. (0->55)
t25 = sin(pkin(3));
t36 = sin(qJ(3,1));
t45 = cos(qJ(3,1));
t46 = cos(qJ(2,1));
t51 = t45 * t46;
t37 = sin(qJ(2,1));
t15 = t37 * t45 * pkin(2) - t46 * pkin(6);
t26 = cos(pkin(3));
t57 = t15 * t26;
t60 = pkin(6) * t37;
t3 = 0.1e1 / (pkin(1) * t57 + (pkin(5) * t60 + (-pkin(1) * t36 + pkin(5) * t51) * pkin(2)) * t25);
t33 = sin(qJ(3,2));
t42 = cos(qJ(3,2));
t43 = cos(qJ(2,2));
t52 = t42 * t43;
t34 = sin(qJ(2,2));
t14 = t34 * t42 * pkin(2) - t43 * pkin(6);
t58 = t14 * t26;
t61 = pkin(6) * t34;
t2 = 0.1e1 / (pkin(1) * t58 + (pkin(5) * t61 + (-pkin(1) * t33 + pkin(5) * t52) * pkin(2)) * t25);
t30 = sin(qJ(3,3));
t39 = cos(qJ(3,3));
t40 = cos(qJ(2,3));
t53 = t39 * t40;
t31 = sin(qJ(2,3));
t13 = t31 * t39 * pkin(2) - t40 * pkin(6);
t59 = t13 * t26;
t62 = pkin(6) * t31;
t1 = 0.1e1 / (pkin(1) * t59 + (pkin(5) * t62 + (-pkin(1) * t30 + pkin(5) * t53) * pkin(2)) * t25);
t63 = pkin(2) * t26;
t56 = t25 * t30;
t55 = t25 * t33;
t54 = t25 * t36;
t47 = cos(qJ(1,1));
t44 = cos(qJ(1,2));
t41 = cos(qJ(1,3));
t38 = sin(qJ(1,1));
t35 = sin(qJ(1,2));
t32 = sin(qJ(1,3));
t29 = legFrame(1,3);
t28 = legFrame(2,3);
t27 = legFrame(3,3);
t21 = cos(t29);
t20 = cos(t28);
t19 = cos(t27);
t18 = sin(t29);
t17 = sin(t28);
t16 = sin(t27);
t12 = t18 * t47 + t38 * t21;
t11 = t17 * t44 + t35 * t20;
t10 = t16 * t41 + t32 * t19;
t9 = t18 * t38 - t47 * t21;
t8 = t17 * t35 - t44 * t20;
t7 = t16 * t32 - t41 * t19;
t4 = [-(t9 * t60 + t12 * t57 + (-t12 * t54 + t9 * t51) * pkin(2)) * t3, -(-t12 * t60 + t9 * t57 + (-t12 * t51 - t9 * t54) * pkin(2)) * t3, (t15 * t25 + t36 * t63) * t3; -(t8 * t61 + t11 * t58 + (-t11 * t55 + t8 * t52) * pkin(2)) * t2, -(-t11 * t61 + t8 * t58 + (-t11 * t52 - t8 * t55) * pkin(2)) * t2, (t14 * t25 + t33 * t63) * t2; -(t7 * t62 + t10 * t59 + (-t10 * t56 + t7 * t53) * pkin(2)) * t1, -(-t10 * t62 + t7 * t59 + (-t10 * t53 - t7 * t56) * pkin(2)) * t1, (t13 * t25 + t30 * t63) * t1;];
Jinv  = t4;
