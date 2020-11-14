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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d3,d4]';
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
% Datum: 2020-08-07 03:40
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRRRR7V1G1P1A3_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(5,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V1G1P1A3_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V1G1P1A3_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRRRR7V1G1P1A3_Jinv: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V1G1P1A3_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V1G1P1A3_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:40:30
% EndTime: 2020-08-07 03:40:30
% DurationCPUTime: 0.09s
% Computational Cost: add. (102->46), mult. (48->41), div. (27->5), fcn. (45->39), ass. (0->27)
t22 = (pkin(4) + pkin(5));
t31 = -2 * t22;
t30 = 2 * t22;
t29 = 0.1e1 / pkin(1) / pkin(2);
t13 = legFrame(3,3) + qJ(1,3);
t14 = legFrame(2,3) + qJ(1,2);
t15 = legFrame(1,3) + qJ(1,1);
t28 = t29 / 0.2e1;
t8 = -qJ(2,3) + t13;
t7 = qJ(2,3) + t13;
t10 = -qJ(2,2) + t14;
t9 = qJ(2,2) + t14;
t12 = -qJ(2,1) + t15;
t11 = qJ(2,1) + t15;
t16 = 0.1e1 / sin(qJ(3,3));
t27 = t16 * t28;
t17 = 0.1e1 / sin(qJ(3,2));
t26 = t17 * t28;
t18 = 0.1e1 / sin(qJ(3,1));
t25 = t18 * t28;
t6 = -qJ(3,1) + t12;
t5 = qJ(3,1) + t11;
t4 = -qJ(3,2) + t10;
t3 = qJ(3,2) + t9;
t2 = -qJ(3,3) + t8;
t1 = qJ(3,3) + t7;
t19 = [(sin(t15) * t31 + (-cos(t6) - cos(t5)) * pkin(2) + (-cos(t12) - cos(t11)) * pkin(1)) * t25, (cos(t15) * t30 + (-sin(t6) - sin(t5)) * pkin(2) + (-sin(t12) - sin(t11)) * pkin(1)) * t25, (-pkin(2) * sin(qJ(2,1) + qJ(3,1)) - sin(qJ(2,1)) * pkin(1)) * t18 * t29; (sin(t14) * t31 + (-cos(t4) - cos(t3)) * pkin(2) + (-cos(t10) - cos(t9)) * pkin(1)) * t26, (cos(t14) * t30 + (-sin(t4) - sin(t3)) * pkin(2) + (-sin(t10) - sin(t9)) * pkin(1)) * t26, (-pkin(2) * sin(qJ(2,2) + qJ(3,2)) - sin(qJ(2,2)) * pkin(1)) * t17 * t29; (sin(t13) * t31 + (-cos(t2) - cos(t1)) * pkin(2) + (-cos(t8) - cos(t7)) * pkin(1)) * t27, (cos(t13) * t30 + (-sin(t2) - sin(t1)) * pkin(2) + (-sin(t8) - sin(t7)) * pkin(1)) * t27, (-pkin(2) * sin(qJ(2,3) + qJ(3,3)) - sin(qJ(2,3)) * pkin(1)) * t16 * t29;];
Jinv  = t19;
