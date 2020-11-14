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
% Datum: 2020-08-06 18:46
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RPRRR6V1G3P3A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(7,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G3P3A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G3P3A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G3P3A1_Jinv: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G3P3A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G3P3A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:45:54
% EndTime: 2020-08-06 18:45:54
% DurationCPUTime: 0.19s
% Computational Cost: add. (111->51), mult. (168->98), div. (18->6), fcn. (177->20), ass. (0->49)
t32 = sin(pkin(7));
t49 = -pkin(6) - pkin(5);
t10 = t49 * t32 - pkin(1);
t33 = cos(pkin(7));
t38 = sin(qJ(1,3));
t44 = cos(qJ(1,3));
t62 = t38 * t32;
t70 = t10 * t44 + pkin(2) * t62 - (pkin(2) * t44 - t38 * t49) * t33;
t40 = sin(qJ(1,2));
t46 = cos(qJ(1,2));
t59 = t40 * t32;
t69 = t10 * t46 + pkin(2) * t59 - (pkin(2) * t46 - t40 * t49) * t33;
t42 = sin(qJ(1,1));
t48 = cos(qJ(1,1));
t56 = t42 * t32;
t68 = t10 * t48 + pkin(2) * t56 - (pkin(2) * t48 - t42 * t49) * t33;
t43 = cos(qJ(3,3));
t13 = t43 * pkin(3) + pkin(2);
t16 = t33 * pkin(1);
t37 = sin(qJ(3,3));
t67 = 0.1e1 / t37 / (t16 + t13);
t45 = cos(qJ(3,2));
t14 = t45 * pkin(3) + pkin(2);
t39 = sin(qJ(3,2));
t66 = 0.1e1 / t39 / (t16 + t14);
t47 = cos(qJ(3,1));
t15 = t47 * pkin(3) + pkin(2);
t41 = sin(qJ(3,1));
t65 = 0.1e1 / t41 / (t16 + t15);
t34 = legFrame(3,2);
t17 = sin(t34);
t64 = t37 * t17;
t20 = cos(t34);
t63 = t37 * t20;
t35 = legFrame(2,2);
t18 = sin(t35);
t61 = t39 * t18;
t21 = cos(t35);
t60 = t39 * t21;
t36 = legFrame(1,2);
t19 = sin(t36);
t58 = t41 * t19;
t22 = cos(t36);
t57 = t41 * t22;
t55 = pkin(3) * (-t44 * t33 + t62) * t43 ^ 2;
t54 = pkin(3) * (-t46 * t33 + t59) * t45 ^ 2;
t53 = pkin(3) * (-t48 * t33 + t56) * t47 ^ 2;
t12 = t16 + pkin(2);
t1 = [(-t22 * t53 + (pkin(3) * t58 - t68 * t22) * t47 + t12 * t58) * t65, (t19 * t53 + (pkin(3) * t57 + t68 * t19) * t47 + t12 * t57) * t65, -t47 * ((t15 * t42 + t48 * t49) * t33 - t10 * t42 + t32 * t48 * t15) * t65; (-t21 * t54 + (pkin(3) * t61 - t69 * t21) * t45 + t12 * t61) * t66, (t18 * t54 + (pkin(3) * t60 + t69 * t18) * t45 + t12 * t60) * t66, -t45 * ((t14 * t40 + t46 * t49) * t33 - t10 * t40 + t32 * t46 * t14) * t66; (-t20 * t55 + (pkin(3) * t64 - t70 * t20) * t43 + t12 * t64) * t67, (t17 * t55 + (pkin(3) * t63 + t70 * t17) * t43 + t12 * t63) * t67, -t43 * ((t13 * t38 + t44 * t49) * t33 - t10 * t38 + t32 * t44 * t13) * t67;];
Jinv  = t1;
