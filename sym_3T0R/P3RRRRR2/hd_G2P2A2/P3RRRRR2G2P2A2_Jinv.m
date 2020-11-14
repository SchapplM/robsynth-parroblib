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

function Jinv = P3RRRRR2G2P2A2_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(2,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G2P2A2_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G2P2A2_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G2P2A2_Jinv: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G2P2A2_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G2P2A2_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:39:24
% EndTime: 2020-08-07 03:39:24
% DurationCPUTime: 0.18s
% Computational Cost: add. (51->42), mult. (102->69), div. (33->8), fcn. (111->36), ass. (0->54)
t37 = cos(qJ(3,1));
t67 = t37 ^ 2;
t34 = cos(qJ(3,2));
t66 = t34 ^ 2;
t31 = cos(qJ(3,3));
t65 = t31 ^ 2;
t64 = -2 * pkin(1);
t24 = sin(qJ(1,3));
t63 = pkin(1) * t24;
t27 = sin(qJ(1,2));
t62 = pkin(1) * t27;
t30 = sin(qJ(1,1));
t61 = pkin(1) * t30;
t22 = sin(qJ(3,3));
t60 = pkin(2) * t22;
t25 = sin(qJ(3,2));
t59 = pkin(2) * t25;
t28 = sin(qJ(3,1));
t58 = pkin(2) * t28;
t57 = 0.1e1 / pkin(2) / pkin(1);
t56 = qJ(2,1) - qJ(3,1);
t55 = qJ(2,1) + qJ(3,1);
t54 = qJ(2,2) - qJ(3,2);
t53 = qJ(2,2) + qJ(3,2);
t52 = qJ(2,3) - qJ(3,3);
t51 = qJ(2,3) + qJ(3,3);
t23 = sin(qJ(2,3));
t32 = cos(qJ(2,3));
t33 = cos(qJ(1,3));
t50 = pkin(2) * (t33 * t23 + t24 * t32) * t65;
t26 = sin(qJ(2,2));
t35 = cos(qJ(2,2));
t36 = cos(qJ(1,2));
t49 = pkin(2) * (t36 * t26 + t27 * t35) * t66;
t29 = sin(qJ(2,1));
t38 = cos(qJ(2,1));
t39 = cos(qJ(1,1));
t48 = pkin(2) * (t39 * t29 + t30 * t38) * t67;
t47 = t32 * t22 * pkin(1);
t46 = t35 * t25 * pkin(1);
t45 = t38 * t28 * pkin(1);
t44 = 0.1e1 / t23 / t65 * t57;
t43 = 0.1e1 / t26 / t66 * t57;
t42 = 0.1e1 / t29 / t67 * t57;
t21 = legFrame(1,2);
t20 = legFrame(2,2);
t19 = legFrame(3,2);
t9 = cos(t21);
t8 = cos(t20);
t7 = cos(t19);
t6 = sin(t21);
t5 = sin(t20);
t4 = sin(t19);
t1 = [(-t9 * t48 + (-t6 * t58 - t9 * t61) * t37 - t6 * t45) * t42, (t6 * t48 + (-t9 * t58 + t6 * t61) * t37 - t9 * t45) * t42, (t39 * t64 + (-cos(qJ(1,1) + t56) - cos(qJ(1,1) + t55)) * pkin(2)) / (sin(t55) + sin(t56)) * t57; (-t8 * t49 + (-t5 * t59 - t8 * t62) * t34 - t5 * t46) * t43, (t5 * t49 + (t5 * t62 - t8 * t59) * t34 - t8 * t46) * t43, (t36 * t64 + (-cos(qJ(1,2) + t54) - cos(qJ(1,2) + t53)) * pkin(2)) / (sin(t53) + sin(t54)) * t57; (-t7 * t50 + (-t4 * t60 - t7 * t63) * t31 - t4 * t47) * t44, (t4 * t50 + (t4 * t63 - t7 * t60) * t31 - t7 * t47) * t44, (t33 * t64 + (-cos(qJ(1,3) + t52) - cos(qJ(1,3) + t51)) * pkin(2)) / (sin(t51) + sin(t52)) * t57;];
Jinv  = t1;
