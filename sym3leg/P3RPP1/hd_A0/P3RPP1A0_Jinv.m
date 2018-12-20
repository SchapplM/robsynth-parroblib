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
% Datum: 2018-12-20 17:50
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Jinv = P3RPP1A0_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPP1A0_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPP1A0_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RPP1A0_Jinv: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPP1A0_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPP1A0_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:50:17
% EndTime: 2018-12-20 17:50:17
% DurationCPUTime: 0.10s
% Computational Cost: add. (135->46), mult. (171->83), div. (9->3), fcn. (78->14), ass. (0->51)
t51 = 2 * pkin(1);
t50 = pkin(1) ^ 2 + 1;
t48 = koppelP(1,1);
t47 = koppelP(2,1);
t46 = koppelP(3,1);
t45 = koppelP(1,2);
t44 = koppelP(2,2);
t43 = koppelP(3,2);
t42 = xP(3);
t41 = cos(qJ(1,1));
t40 = cos(qJ(1,2));
t39 = cos(qJ(1,3));
t38 = sin(qJ(1,1));
t37 = sin(qJ(1,2));
t36 = sin(qJ(1,3));
t35 = pkin(1) + qJ(3,1);
t34 = pkin(1) + qJ(3,2);
t33 = pkin(1) + qJ(3,3);
t32 = legFrame(1,3);
t31 = legFrame(2,3);
t30 = legFrame(3,3);
t29 = cos(t42);
t28 = sin(t42);
t27 = cos(t32);
t26 = cos(t31);
t25 = cos(t30);
t24 = sin(t32);
t23 = sin(t31);
t22 = sin(t30);
t21 = qJ(2,1) * t48 + t35 * t45;
t20 = qJ(2,2) * t47 + t34 * t44;
t19 = qJ(2,3) * t46 + t33 * t43;
t18 = -qJ(2,1) * t45 + t35 * t48;
t17 = -qJ(2,2) * t44 + t34 * t47;
t16 = -qJ(2,3) * t43 + t33 * t46;
t15 = 1 / (qJ(2,1) ^ 2 + (t51 + qJ(3,1)) * qJ(3,1) + t50);
t14 = 1 / (qJ(2,2) ^ 2 + (t51 + qJ(3,2)) * qJ(3,2) + t50);
t13 = 1 / (qJ(2,3) ^ 2 + (t51 + qJ(3,3)) * qJ(3,3) + t50);
t12 = t38 * qJ(2,1) + t35 * t41;
t11 = t37 * qJ(2,2) + t34 * t40;
t10 = t36 * qJ(2,3) + t33 * t39;
t9 = t41 * qJ(2,1) - t38 * t35;
t8 = t40 * qJ(2,2) - t37 * t34;
t7 = t39 * qJ(2,3) - t36 * t33;
t6 = t38 * t18 - t21 * t41;
t5 = t37 * t17 - t20 * t40;
t4 = t36 * t16 - t19 * t39;
t3 = t18 * t41 + t21 * t38;
t2 = t17 * t40 + t20 * t37;
t1 = t16 * t39 + t19 * t36;
t49 = [(-t24 * t12 + t9 * t27) * t15 (t12 * t27 + t24 * t9) * t15 ((t6 * t28 + t3 * t29) * t27 - t24 * (-t3 * t28 + t6 * t29)) * t15; (-t23 * t11 + t8 * t26) * t14 (t11 * t26 + t23 * t8) * t14 ((t2 * t29 + t5 * t28) * t26 - t23 * (-t2 * t28 + t5 * t29)) * t14; (-t22 * t10 + t7 * t25) * t13 (t10 * t25 + t22 * t7) * t13 ((t1 * t29 + t4 * t28) * t25 - t22 * (-t1 * t28 + t4 * t29)) * t13;];
Jinv  = t49;
