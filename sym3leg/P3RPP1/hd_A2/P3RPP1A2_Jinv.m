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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
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
% Datum: 2018-12-20 17:53
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Jinv = P3RPP1A2_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPP1A2_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPP1A2_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RPP1A2_Jinv: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPP1A2_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPP1A2_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:53:04
% EndTime: 2018-12-20 17:53:04
% DurationCPUTime: 0.11s
% Computational Cost: add. (171->55), mult. (243->86), div. (9->3), fcn. (78->14), ass. (0->57)
t60 = 2 * pkin(1);
t59 = pkin(1) ^ 2 + 1;
t58 = qJ(2,1) * (pkin(1) + qJ(3,1));
t57 = qJ(2,2) * (pkin(1) + qJ(3,2));
t56 = qJ(2,3) * (pkin(1) + qJ(3,3));
t54 = koppelP(1,1);
t53 = koppelP(2,1);
t52 = koppelP(3,1);
t51 = koppelP(1,2);
t50 = koppelP(2,2);
t49 = koppelP(3,2);
t48 = qJ(2,1) ^ 2;
t47 = qJ(2,2) ^ 2;
t46 = qJ(2,3) ^ 2;
t45 = xP(3);
t44 = cos(qJ(1,1));
t43 = cos(qJ(1,2));
t42 = cos(qJ(1,3));
t41 = sin(qJ(1,1));
t40 = sin(qJ(1,2));
t39 = sin(qJ(1,3));
t35 = legFrame(1,3);
t34 = legFrame(2,3);
t33 = legFrame(3,3);
t32 = 1 + t48;
t31 = 1 + t47;
t30 = 1 + t46;
t29 = cos(t45);
t28 = sin(t45);
t27 = cos(t35);
t26 = cos(t34);
t25 = cos(t33);
t24 = sin(t35);
t23 = sin(t34);
t22 = sin(t33);
t21 = 1 / (t48 + (t60 + qJ(3,1)) * qJ(3,1) + t59);
t20 = 1 / (t47 + (t60 + qJ(3,2)) * qJ(3,2) + t59);
t19 = 1 / (t46 + (t60 + qJ(3,3)) * qJ(3,3) + t59);
t18 = t48 * t54 + t51 * t58 + t54;
t17 = -t48 * t51 + t54 * t58 - t51;
t16 = t47 * t53 + t50 * t57 + t53;
t15 = -t47 * t50 + t53 * t57 - t50;
t14 = t46 * t52 + t49 * t56 + t52;
t13 = -t46 * t49 + t52 * t56 - t49;
t12 = -t32 * t44 + t41 * t58;
t11 = -t31 * t43 + t40 * t57;
t10 = -t30 * t42 + t39 * t56;
t9 = t41 * t32 + t44 * t58;
t8 = t40 * t31 + t43 * t57;
t7 = t39 * t30 + t42 * t56;
t6 = t17 * t41 - t18 * t44;
t5 = t17 * t44 + t18 * t41;
t4 = t15 * t40 - t16 * t43;
t3 = t15 * t43 + t16 * t40;
t2 = t13 * t39 - t14 * t42;
t1 = t13 * t42 + t14 * t39;
t36 = [(-t12 * t24 + t9 * t27) * t21 (t12 * t27 + t9 * t24) * t21 ((-t28 * t5 + t6 * t29) * t27 + (t28 * t6 + t5 * t29) * t24) * t21; (-t11 * t23 + t8 * t26) * t20 (t11 * t26 + t8 * t23) * t20 ((-t28 * t3 + t4 * t29) * t26 + (t28 * t4 + t3 * t29) * t23) * t20; (-t10 * t22 + t7 * t25) * t19 (t10 * t25 + t7 * t22) * t19 ((-t28 * t1 + t2 * t29) * t25 + (t1 * t29 + t28 * t2) * t22) * t19;];
Jinv  = t36;
