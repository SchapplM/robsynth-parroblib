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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2020-08-07 03:39
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRRRR2G3P3A2_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(2,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G3P3A2_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G3P3A2_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G3P3A2_Jinv: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G3P3A2_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G3P3A2_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:39:42
% EndTime: 2020-08-07 03:39:42
% DurationCPUTime: 0.18s
% Computational Cost: add. (30->27), mult. (105->77), div. (36->11), fcn. (120->24), ass. (0->50)
t37 = cos(qJ(3,1));
t63 = t37 ^ 2;
t34 = cos(qJ(3,2));
t62 = t34 ^ 2;
t31 = cos(qJ(3,3));
t61 = t31 ^ 2;
t33 = cos(qJ(1,3));
t60 = pkin(1) * t33;
t36 = cos(qJ(1,2));
t59 = pkin(1) * t36;
t39 = cos(qJ(1,1));
t58 = pkin(1) * t39;
t22 = sin(qJ(3,3));
t57 = pkin(2) * t22;
t25 = sin(qJ(3,2));
t56 = pkin(2) * t25;
t28 = sin(qJ(3,1));
t55 = pkin(2) * t28;
t54 = 0.1e1 / pkin(2) / pkin(1);
t23 = sin(qJ(2,3));
t24 = sin(qJ(1,3));
t32 = cos(qJ(2,3));
t53 = pkin(2) * (t24 * t23 - t33 * t32) * t61;
t26 = sin(qJ(2,2));
t27 = sin(qJ(1,2));
t35 = cos(qJ(2,2));
t52 = pkin(2) * (t27 * t26 - t36 * t35) * t62;
t29 = sin(qJ(2,1));
t30 = sin(qJ(1,1));
t38 = cos(qJ(2,1));
t51 = pkin(2) * (t30 * t29 - t39 * t38) * t63;
t50 = pkin(1) * t22 * t32;
t49 = pkin(1) * t25 * t35;
t48 = pkin(1) * t28 * t38;
t47 = 0.1e1 / t23 * t54;
t46 = 0.1e1 / t26 * t54;
t45 = 0.1e1 / t29 * t54;
t44 = 0.1e1 / t61 * t47;
t43 = 0.1e1 / t62 * t46;
t42 = 0.1e1 / t63 * t45;
t21 = legFrame(1,2);
t20 = legFrame(2,2);
t19 = legFrame(3,2);
t9 = cos(t21);
t8 = cos(t20);
t7 = cos(t19);
t6 = sin(t21);
t5 = sin(t20);
t4 = sin(t19);
t1 = [(t9 * t51 + (-t6 * t55 - t9 * t58) * t37 - t6 * t48) * t42, (-t6 * t51 + (-t9 * t55 + t6 * t58) * t37 - t9 * t48) * t42, (pkin(2) * (t39 * t29 + t30 * t38) * t37 + t30 * pkin(1)) / t37 * t45; (t8 * t52 + (-t5 * t56 - t8 * t59) * t34 - t5 * t49) * t43, (-t5 * t52 + (t5 * t59 - t8 * t56) * t34 - t8 * t49) * t43, (pkin(2) * (t36 * t26 + t27 * t35) * t34 + t27 * pkin(1)) / t34 * t46; (t7 * t53 + (-t4 * t57 - t7 * t60) * t31 - t4 * t50) * t44, (-t4 * t53 + (t4 * t60 - t7 * t57) * t31 - t7 * t50) * t44, (pkin(2) * (t33 * t23 + t24 * t32) * t31 + t24 * pkin(1)) / t31 * t47;];
Jinv  = t1;
