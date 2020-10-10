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
%   pkin=[a3,a4,d4]';
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
% Datum: 2020-08-06 19:58
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRPRR1V1G2P2A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR1V1G2P2A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR1V1G2P2A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RRPRR1V1G2P2A1_Jinv: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR1V1G2P2A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR1V1G2P2A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:58:16
% EndTime: 2020-08-06 19:58:16
% DurationCPUTime: 0.08s
% Computational Cost: add. (48->16), mult. (45->39), div. (9->3), fcn. (45->18), ass. (0->32)
t28 = pkin(2) + pkin(1);
t37 = sin(qJ(2,3)) * t28;
t36 = sin(qJ(2,2)) * t28;
t35 = sin(qJ(2,1)) * t28;
t34 = cos(qJ(2,3)) * t28;
t33 = cos(qJ(2,2)) * t28;
t32 = cos(qJ(2,1)) * t28;
t10 = pkin(3) + qJ(3,3);
t17 = sin(qJ(1,3));
t23 = cos(qJ(1,3));
t31 = -t10 * t23 + t17 * t34;
t11 = pkin(3) + qJ(3,2);
t19 = sin(qJ(1,2));
t25 = cos(qJ(1,2));
t30 = -t11 * t25 + t19 * t33;
t12 = pkin(3) + qJ(3,1);
t21 = sin(qJ(1,1));
t27 = cos(qJ(1,1));
t29 = -t12 * t27 + t21 * t32;
t15 = legFrame(1,2);
t14 = legFrame(2,2);
t13 = legFrame(3,2);
t9 = 0.1e1 / t12;
t8 = 0.1e1 / t11;
t7 = 0.1e1 / t10;
t6 = cos(t15);
t5 = cos(t14);
t4 = cos(t13);
t3 = sin(t15);
t2 = sin(t14);
t1 = sin(t13);
t16 = [(t29 * t6 + t3 * t35) * t9, (-t29 * t3 + t6 * t35) * t9, (t21 * t12 + t27 * t32) * t9; (t2 * t36 + t30 * t5) * t8, (-t30 * t2 + t5 * t36) * t8, (t19 * t11 + t25 * t33) * t8; (t1 * t37 + t31 * t4) * t7, (-t31 * t1 + t4 * t37) * t7, (t17 * t10 + t23 * t34) * t7;];
Jinv  = t16;
