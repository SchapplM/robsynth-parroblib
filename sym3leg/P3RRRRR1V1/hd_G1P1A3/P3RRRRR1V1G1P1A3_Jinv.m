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

function Jinv = P3RRRRR1V1G1P1A3_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(4,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR1V1G1P1A3_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR1V1G1P1A3_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRRRR1V1G1P1A3_Jinv: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR1V1G1P1A3_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR1V1G1P1A3_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:33:07
% EndTime: 2020-08-07 03:33:07
% DurationCPUTime: 0.07s
% Computational Cost: add. (30->24), mult. (48->41), div. (27->5), fcn. (63->27), ass. (0->29)
t36 = 0.1e1 / pkin(2) / pkin(3);
t16 = sin(qJ(3,3));
t35 = 0.1e1 / t16 * t36;
t19 = sin(qJ(3,2));
t34 = 0.1e1 / t19 * t36;
t22 = sin(qJ(3,1));
t33 = 0.1e1 / t22 * t36;
t17 = sin(qJ(2,3));
t32 = ((pkin(3) * cos(qJ(3,3)) + pkin(2)) * cos(qJ(2,3)) - t16 * t17 * pkin(3)) * t35;
t20 = sin(qJ(2,2));
t31 = ((pkin(3) * cos(qJ(3,2)) + pkin(2)) * cos(qJ(2,2)) - t19 * t20 * pkin(3)) * t34;
t23 = sin(qJ(2,1));
t30 = ((pkin(3) * cos(qJ(3,1)) + pkin(2)) * cos(qJ(2,1)) - t22 * t23 * pkin(3)) * t33;
t27 = cos(qJ(1,1));
t26 = cos(qJ(1,2));
t25 = cos(qJ(1,3));
t24 = sin(qJ(1,1));
t21 = sin(qJ(1,2));
t18 = sin(qJ(1,3));
t15 = legFrame(1,3);
t14 = legFrame(2,3);
t13 = legFrame(3,3);
t9 = cos(t15);
t8 = cos(t14);
t7 = cos(t13);
t6 = sin(t15);
t5 = sin(t14);
t4 = sin(t13);
t1 = [-(-t6 * t24 + t9 * t27) * t30, -(t9 * t24 + t6 * t27) * t30, (pkin(3) * sin(qJ(2,1) + qJ(3,1)) + t23 * pkin(2)) * t33; -(-t5 * t21 + t8 * t26) * t31, -(t8 * t21 + t5 * t26) * t31, (pkin(3) * sin(qJ(2,2) + qJ(3,2)) + t20 * pkin(2)) * t34; -(-t4 * t18 + t7 * t25) * t32, -(t7 * t18 + t4 * t25) * t32, (pkin(3) * sin(qJ(2,3) + qJ(3,3)) + t17 * pkin(2)) * t35;];
Jinv  = t1;
