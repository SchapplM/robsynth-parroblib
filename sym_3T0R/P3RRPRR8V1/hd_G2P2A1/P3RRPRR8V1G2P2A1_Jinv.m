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
%   pkin=[a3,a4,d1,d4,theta3]';
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
% Datum: 2020-08-06 21:01
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRPRR8V1G2P2A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(5,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V1G2P2A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V1G2P2A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRPRR8V1G2P2A1_Jinv: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V1G2P2A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V1G2P2A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 21:01:38
% EndTime: 2020-08-06 21:01:39
% DurationCPUTime: 0.08s
% Computational Cost: add. (63->28), mult. (93->53), div. (9->3), fcn. (93->23), ass. (0->36)
t36 = pkin(2) * sin(pkin(5));
t35 = cos(qJ(1,1));
t34 = cos(qJ(2,1));
t33 = cos(qJ(1,2));
t32 = cos(qJ(2,2));
t31 = cos(qJ(1,3));
t30 = cos(qJ(2,3));
t29 = sin(qJ(1,1));
t28 = sin(qJ(2,1));
t27 = sin(qJ(1,2));
t26 = sin(qJ(2,2));
t25 = sin(qJ(1,3));
t24 = sin(qJ(2,3));
t23 = legFrame(1,2);
t22 = legFrame(2,2);
t21 = legFrame(3,2);
t20 = pkin(4) + qJ(3,1);
t19 = pkin(4) + qJ(3,2);
t18 = pkin(4) + qJ(3,3);
t16 = 0.1e1 / t20;
t15 = 0.1e1 / t19;
t14 = 0.1e1 / t18;
t13 = cos(t23);
t12 = cos(t22);
t11 = cos(t21);
t10 = sin(t23);
t9 = sin(t22);
t8 = sin(t21);
t7 = cos(pkin(5)) * pkin(2) + pkin(1);
t6 = t28 * t7 + t34 * t36;
t5 = t26 * t7 + t32 * t36;
t4 = t24 * t7 + t30 * t36;
t3 = -t35 * t20 + (-t28 * t36 + t34 * t7) * t29;
t2 = -t33 * t19 + (-t26 * t36 + t32 * t7) * t27;
t1 = -t31 * t18 + (-t24 * t36 + t30 * t7) * t25;
t17 = [(t10 * t6 + t3 * t13) * t16, (-t3 * t10 + t13 * t6) * t16, (t29 * t20 + (pkin(2) * cos(qJ(2,1) + pkin(5)) + pkin(1) * t34) * t35) * t16; (t2 * t12 + t9 * t5) * t15, (t12 * t5 - t2 * t9) * t15, (t27 * t19 + (pkin(2) * cos(qJ(2,2) + pkin(5)) + pkin(1) * t32) * t33) * t15; (t1 * t11 + t8 * t4) * t14, (-t1 * t8 + t11 * t4) * t14, (t25 * t18 + (pkin(2) * cos(qJ(2,3) + pkin(5)) + pkin(1) * t30) * t31) * t14;];
Jinv  = t17;
