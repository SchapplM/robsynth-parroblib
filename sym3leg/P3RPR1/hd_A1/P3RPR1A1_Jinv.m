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
% Datum: 2019-05-03 14:59
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RPR1A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(4,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RPR1A1_Jinv: qJ has to be [2x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPR1A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RPR1A1_Jinv: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPR1A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPR1A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:59:23
% EndTime: 2019-05-03 14:59:23
% DurationCPUTime: 0.07s
% Computational Cost: add. (27->21), mult. (54->51), div. (9->3), fcn. (66->14), ass. (0->34)
t33 = koppelP(1,1);
t32 = koppelP(2,1);
t31 = koppelP(3,1);
t30 = koppelP(1,2);
t29 = koppelP(2,2);
t28 = koppelP(3,2);
t27 = 0.1e1 / qJ(2,1);
t26 = 0.1e1 / qJ(2,2);
t25 = 0.1e1 / qJ(2,3);
t24 = xP(3);
t23 = cos(qJ(1,1));
t22 = cos(qJ(1,2));
t21 = cos(qJ(1,3));
t20 = sin(qJ(1,1));
t19 = sin(qJ(1,2));
t18 = sin(qJ(1,3));
t17 = legFrame(1,3);
t16 = legFrame(2,3);
t15 = legFrame(3,3);
t14 = cos(t24);
t13 = sin(t24);
t12 = cos(t17);
t11 = cos(t16);
t10 = cos(t15);
t9 = sin(t17);
t8 = sin(t16);
t7 = sin(t15);
t6 = t20 * t30 + t23 * t33;
t5 = t20 * t33 - t23 * t30;
t4 = t19 * t29 + t22 * t32;
t3 = t19 * t32 - t22 * t29;
t2 = t18 * t28 + t21 * t31;
t1 = t18 * t31 - t21 * t28;
t34 = [t27 * (t12 * t23 - t9 * t20), t27 * (t12 * t20 + t9 * t23), ((-t13 * t6 + t5 * t14) * t12 + t9 * (t13 * t5 + t6 * t14)) * t27; t26 * (t11 * t22 - t8 * t19), t26 * (t11 * t19 + t8 * t22), ((-t13 * t4 + t3 * t14) * t11 + t8 * (t13 * t3 + t4 * t14)) * t26; t25 * (t10 * t21 - t7 * t18), t25 * (t10 * t18 + t7 * t21), ((t1 * t14 - t13 * t2) * t10 + t7 * (t13 * t1 + t2 * t14)) * t25;];
Jinv  = t34;
