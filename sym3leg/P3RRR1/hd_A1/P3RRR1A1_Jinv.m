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
% Datum: 2019-05-03 15:40
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRR1A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(5,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RRR1A1_Jinv: qJ has to be [2x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRR1A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRR1A1_Jinv: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRR1A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRR1A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 15:40:06
% EndTime: 2019-05-03 15:40:06
% DurationCPUTime: 0.13s
% Computational Cost: add. (72->32), mult. (72->57), div. (18->4), fcn. (102->20), ass. (0->38)
t29 = qJ(1,1) + qJ(2,1);
t15 = sin(t29);
t18 = cos(t29);
t40 = 0.1e1 / pkin(1);
t43 = 0.1e1 / (t15 * cos(qJ(1,1)) - sin(qJ(1,1)) * t18) * t40;
t27 = qJ(1,3) + qJ(2,3);
t13 = sin(t27);
t16 = cos(t27);
t42 = 0.1e1 / (t13 * cos(qJ(1,3)) - sin(qJ(1,3)) * t16) * t40;
t28 = qJ(1,2) + qJ(2,2);
t14 = sin(t28);
t17 = cos(t28);
t41 = 0.1e1 / (t14 * cos(qJ(1,2)) - sin(qJ(1,2)) * t17) * t40;
t39 = koppelP(1,1);
t38 = koppelP(2,1);
t37 = koppelP(3,1);
t36 = koppelP(1,2);
t35 = koppelP(2,2);
t34 = koppelP(3,2);
t33 = xP(3);
t32 = legFrame(1,3);
t31 = legFrame(2,3);
t30 = legFrame(3,3);
t26 = cos(t33);
t25 = sin(t33);
t24 = cos(t32);
t23 = cos(t31);
t22 = cos(t30);
t21 = sin(t32);
t20 = sin(t31);
t19 = sin(t30);
t12 = -t25 * t36 + t26 * t39;
t11 = -t25 * t35 + t26 * t38;
t10 = -t25 * t34 + t26 * t37;
t9 = t25 * t39 + t26 * t36;
t8 = t25 * t38 + t26 * t35;
t7 = t25 * t37 + t26 * t34;
t1 = [(-t15 * t21 + t24 * t18) * t43, (t24 * t15 + t21 * t18) * t43, -((-t21 * t12 + t9 * t24) * t18 - t15 * (t12 * t24 + t21 * t9)) * t43; -(t14 * t20 - t23 * t17) * t41, (t23 * t14 + t20 * t17) * t41, -((-t20 * t11 + t8 * t23) * t17 - t14 * (t11 * t23 + t20 * t8)) * t41; -(t13 * t19 - t22 * t16) * t42, (t22 * t13 + t19 * t16) * t42, -((-t19 * t10 + t7 * t22) * t16 - t13 * (t10 * t22 + t19 * t7)) * t42;];
Jinv  = t1;
