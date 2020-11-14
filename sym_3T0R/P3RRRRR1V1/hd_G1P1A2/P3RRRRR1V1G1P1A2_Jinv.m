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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1]';
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
% Datum: 2020-08-07 03:33
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRRRR1V1G1P1A2_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(4,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR1V1G1P1A2_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR1V1G1P1A2_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRRRR1V1G1P1A2_Jinv: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR1V1G1P1A2_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR1V1G1P1A2_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:32:56
% EndTime: 2020-08-07 03:32:57
% DurationCPUTime: 0.11s
% Computational Cost: add. (48->21), mult. (6->12), div. (18->4), fcn. (24->18), ass. (0->20)
t13 = 0.1e1 / pkin(2);
t25 = t13 / sin(qJ(3,3));
t24 = t13 / sin(qJ(3,2));
t23 = t13 / sin(qJ(3,1));
t22 = legFrame(1,3) + qJ(1,1);
t21 = qJ(1,2) + legFrame(2,3);
t20 = qJ(1,3) + legFrame(3,3);
t19 = qJ(2,1) + qJ(3,1);
t18 = qJ(2,2) + qJ(3,2);
t17 = qJ(2,3) + qJ(3,3);
t16 = t25 / 0.2e1;
t15 = t24 / 0.2e1;
t14 = t23 / 0.2e1;
t6 = -t19 + t22;
t5 = t19 + t22;
t4 = -t18 + t21;
t3 = t18 + t21;
t2 = -t17 + t20;
t1 = t17 + t20;
t7 = [(cos(t6) + cos(t5)) * t14, (sin(t5) + sin(t6)) * t14, -sin(t19) * t23; (cos(t4) + cos(t3)) * t15, (sin(t3) + sin(t4)) * t15, -sin(t18) * t24; (cos(t2) + cos(t1)) * t16, (sin(t1) + sin(t2)) * t16, -sin(t17) * t25;];
Jinv  = t7;
