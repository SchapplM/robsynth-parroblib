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
%   pkin=[a3,a4,d1,d4]';
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
% Datum: 2020-08-06 19:05
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRPRR12V1G1P1A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(4,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G1P1A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G1P1A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G1P1A1_Jinv: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G1P1A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G1P1A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:05:10
% EndTime: 2020-08-06 19:05:10
% DurationCPUTime: 0.06s
% Computational Cost: add. (48->19), mult. (66->45), div. (9->3), fcn. (66->18), ass. (0->35)
t34 = 0.1e1 / qJ(3,1);
t33 = 0.1e1 / qJ(3,2);
t32 = 0.1e1 / qJ(3,3);
t31 = pkin(1) + pkin(2);
t30 = cos(qJ(1,1));
t29 = cos(qJ(2,1));
t28 = cos(qJ(1,2));
t27 = cos(qJ(2,2));
t26 = cos(qJ(1,3));
t25 = cos(qJ(2,3));
t24 = sin(qJ(1,1));
t23 = sin(qJ(2,1));
t22 = sin(qJ(1,2));
t21 = sin(qJ(2,2));
t20 = sin(qJ(1,3));
t19 = sin(qJ(2,3));
t18 = legFrame(1,3);
t17 = legFrame(2,3);
t16 = legFrame(3,3);
t15 = cos(t18);
t14 = cos(t17);
t13 = cos(t16);
t12 = sin(t18);
t11 = sin(t17);
t10 = sin(t16);
t9 = t23 * qJ(3,1) + t31 * t29;
t8 = t21 * qJ(3,2) + t31 * t27;
t7 = t19 * qJ(3,3) + t31 * t25;
t6 = -t24 * pkin(4) + t9 * t30;
t5 = -t22 * pkin(4) + t8 * t28;
t4 = -t20 * pkin(4) + t7 * t26;
t3 = t30 * pkin(4) + t9 * t24;
t2 = t28 * pkin(4) + t8 * t22;
t1 = t26 * pkin(4) + t7 * t20;
t35 = [(-t3 * t12 + t15 * t6) * t34, (t12 * t6 + t3 * t15) * t34, (-t29 * qJ(3,1) + t31 * t23) * t34; (-t2 * t11 + t14 * t5) * t33, (t11 * t5 + t2 * t14) * t33, (-t27 * qJ(3,2) + t31 * t21) * t33; (-t1 * t10 + t13 * t4) * t32, (t1 * t13 + t10 * t4) * t32, (-t25 * qJ(3,3) + t31 * t19) * t32;];
Jinv  = t35;
