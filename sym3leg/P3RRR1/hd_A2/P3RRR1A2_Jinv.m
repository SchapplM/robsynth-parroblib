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
% qJ [2x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
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
% Datum: 2019-05-03 15:41
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRR1G1P1A2_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(5,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RRR1G1P1A2_Jinv: qJ has to be [2x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRR1G1P1A2_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRR1G1P1A2_Jinv: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRR1G1P1A2_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRR1G1P1A2_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 15:41:38
% EndTime: 2019-05-03 15:41:38
% DurationCPUTime: 0.16s
% Computational Cost: add. (108->57), mult. (153->116), div. (27->5), fcn. (168->20), ass. (0->50)
t54 = 0.1e1 / pkin(2) / pkin(1);
t30 = qJ(1,3) + qJ(2,3);
t16 = sin(t30);
t19 = cos(t30);
t36 = sin(qJ(1,3));
t39 = cos(qJ(1,3));
t53 = 0.1e1 / (-t16 * t39 + t36 * t19) * t54;
t31 = qJ(1,2) + qJ(2,2);
t17 = sin(t31);
t20 = cos(t31);
t37 = sin(qJ(1,2));
t40 = cos(qJ(1,2));
t52 = 0.1e1 / (-t17 * t40 + t37 * t20) * t54;
t32 = qJ(1,1) + qJ(2,1);
t18 = sin(t32);
t21 = cos(t32);
t38 = sin(qJ(1,1));
t41 = cos(qJ(1,1));
t51 = 0.1e1 / (-t18 * t41 + t38 * t21) * t54;
t48 = koppelP(1,1);
t47 = koppelP(2,1);
t46 = koppelP(3,1);
t45 = koppelP(1,2);
t44 = koppelP(2,2);
t43 = koppelP(3,2);
t42 = xP(3);
t35 = legFrame(1,3);
t34 = legFrame(2,3);
t33 = legFrame(3,3);
t29 = cos(t42);
t28 = sin(t42);
t27 = cos(t35);
t26 = cos(t34);
t25 = cos(t33);
t24 = sin(t35);
t23 = sin(t34);
t22 = sin(t33);
t15 = t38 * t45 + t41 * t48;
t14 = t38 * t48 - t41 * t45;
t13 = t37 * t44 + t40 * t47;
t12 = t37 * t47 - t40 * t44;
t11 = t36 * t43 + t39 * t46;
t10 = t36 * t46 - t39 * t43;
t9 = -t28 * t45 + t29 * t48;
t8 = -t28 * t44 + t29 * t47;
t7 = -t28 * t43 + t29 * t46;
t6 = t28 * t48 + t29 * t45;
t5 = t28 * t47 + t29 * t44;
t4 = t28 * t46 + t29 * t43;
t1 = [(pkin(1) * (-t38 * t24 + t27 * t41) + (-t24 * t18 + t27 * t21) * pkin(2)) * t51, (pkin(1) * (t24 * t41 + t27 * t38) + (t27 * t18 + t24 * t21) * pkin(2)) * t51, (((t14 * t29 - t28 * t15) * t27 + t24 * (t28 * t14 + t15 * t29)) * pkin(1) + (-(-t24 * t9 + t6 * t27) * t21 + (t24 * t6 + t9 * t27) * t18) * pkin(2)) * t51; (pkin(1) * (-t37 * t23 + t26 * t40) + (-t23 * t17 + t26 * t20) * pkin(2)) * t52, (pkin(1) * (t23 * t40 + t26 * t37) + (t26 * t17 + t23 * t20) * pkin(2)) * t52, (((t12 * t29 - t28 * t13) * t26 + t23 * (t28 * t12 + t13 * t29)) * pkin(1) + (-(-t23 * t8 + t5 * t26) * t20 + (t23 * t5 + t8 * t26) * t17) * pkin(2)) * t52; (pkin(1) * (-t36 * t22 + t25 * t39) + (-t22 * t16 + t25 * t19) * pkin(2)) * t53, (pkin(1) * (t22 * t39 + t25 * t36) + (t25 * t16 + t22 * t19) * pkin(2)) * t53, (((t10 * t29 - t28 * t11) * t25 + t22 * (t28 * t10 + t11 * t29)) * pkin(1) + (-(-t22 * t7 + t4 * t25) * t19 + (t22 * t4 + t7 * t25) * t16) * pkin(2)) * t53;];
Jinv  = t1;
