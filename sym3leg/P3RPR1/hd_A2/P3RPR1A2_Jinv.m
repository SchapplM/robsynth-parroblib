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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
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
% Datum: 2019-05-03 15:00
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RPR1A2_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(4,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RPR1A2_Jinv: qJ has to be [2x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPR1A2_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RPR1A2_Jinv: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPR1A2_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPR1A2_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 15:00:53
% EndTime: 2019-05-03 15:00:53
% DurationCPUTime: 0.09s
% Computational Cost: add. (99->34), mult. (126->75), div. (9->3), fcn. (78->14), ass. (0->47)
t46 = koppelP(1,1);
t45 = koppelP(2,1);
t44 = koppelP(3,1);
t43 = koppelP(1,2);
t42 = koppelP(2,2);
t41 = koppelP(3,2);
t40 = 0.1e1 / qJ(2,1);
t39 = 0.1e1 / qJ(2,2);
t38 = 0.1e1 / qJ(2,3);
t37 = xP(3);
t36 = pkin(1) + pkin(2);
t35 = cos(qJ(1,1));
t34 = cos(qJ(1,2));
t33 = cos(qJ(1,3));
t32 = sin(qJ(1,1));
t31 = sin(qJ(1,2));
t30 = sin(qJ(1,3));
t29 = legFrame(1,3);
t28 = legFrame(2,3);
t27 = legFrame(3,3);
t26 = cos(t37);
t25 = sin(t37);
t24 = cos(t29);
t23 = cos(t28);
t22 = cos(t27);
t21 = sin(t29);
t20 = sin(t28);
t19 = sin(t27);
t18 = qJ(2,1) * t46 + t36 * t43;
t17 = qJ(2,2) * t45 + t36 * t42;
t16 = qJ(2,3) * t44 + t36 * t41;
t15 = -qJ(2,1) * t43 + t36 * t46;
t14 = -qJ(2,2) * t42 + t36 * t45;
t13 = -qJ(2,3) * t41 + t36 * t44;
t12 = t32 * qJ(2,1) + t36 * t35;
t11 = t31 * qJ(2,2) + t36 * t34;
t10 = t30 * qJ(2,3) + t36 * t33;
t9 = -t35 * qJ(2,1) + t32 * t36;
t8 = -t34 * qJ(2,2) + t31 * t36;
t7 = -t33 * qJ(2,3) + t30 * t36;
t6 = t15 * t32 - t18 * t35;
t5 = t14 * t31 - t17 * t34;
t4 = t13 * t30 - t16 * t33;
t3 = t15 * t35 + t32 * t18;
t2 = t14 * t34 + t31 * t17;
t1 = t13 * t33 + t30 * t16;
t47 = [(t12 * t24 - t9 * t21) * t40, (t21 * t12 + t9 * t24) * t40, ((-t25 * t3 + t6 * t26) * t24 + t21 * (t25 * t6 + t3 * t26)) * t40; (t11 * t23 - t8 * t20) * t39, (t20 * t11 + t8 * t23) * t39, ((-t25 * t2 + t5 * t26) * t23 + t20 * (t2 * t26 + t25 * t5)) * t39; (t10 * t22 - t7 * t19) * t38, (t19 * t10 + t7 * t22) * t38, ((-t25 * t1 + t4 * t26) * t22 + t19 * (t1 * t26 + t25 * t4)) * t38;];
Jinv  = t47;
