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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2020-08-06 18:24
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RPRRR12V1G1P1A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(6,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G1P1A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G1P1A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G1P1A1_Jinv: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G1P1A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G1P1A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:24:09
% EndTime: 2020-08-06 18:24:09
% DurationCPUTime: 0.11s
% Computational Cost: add. (60->16), mult. (54->33), div. (9->6), fcn. (60->18), ass. (0->32)
t31 = cos(qJ(1,1));
t30 = cos(qJ(1,2));
t29 = cos(qJ(1,3));
t28 = sin(qJ(1,1));
t27 = sin(qJ(3,1));
t26 = sin(qJ(1,2));
t25 = sin(qJ(3,2));
t24 = sin(qJ(1,3));
t23 = sin(qJ(3,3));
t22 = legFrame(1,3);
t21 = legFrame(2,3);
t20 = legFrame(3,3);
t19 = pkin(1) + pkin(5) + pkin(6);
t18 = cos(t22);
t17 = cos(t21);
t16 = cos(t20);
t15 = sin(t22);
t14 = sin(t21);
t13 = sin(t20);
t12 = pkin(3) * t27 + qJ(2,1);
t11 = pkin(3) * t25 + qJ(2,2);
t10 = pkin(3) * t23 + qJ(2,3);
t9 = 0.1e1 / t12;
t8 = 0.1e1 / t11;
t7 = 0.1e1 / t10;
t6 = -t12 * t31 + t19 * t28;
t5 = -t11 * t30 + t19 * t26;
t4 = -t10 * t29 + t19 * t24;
t3 = t12 * t28 + t19 * t31;
t2 = t11 * t26 + t19 * t30;
t1 = t10 * t24 + t19 * t29;
t32 = [(-t15 * t6 + t18 * t3) * t9, (t15 * t3 + t18 * t6) * t9, cos(qJ(3,1)) / t27; (-t14 * t5 + t17 * t2) * t8, (t14 * t2 + t17 * t5) * t8, cos(qJ(3,2)) / t25; (t1 * t16 - t13 * t4) * t7, (t1 * t13 + t16 * t4) * t7, cos(qJ(3,3)) / t23;];
Jinv  = t32;
