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
%   pkin=[a2,a3,d2,d3]';
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
% Datum: 2018-12-20 17:43
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Jinv = P3PRR1A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(4,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PRR1A1_Jinv: qJ has to be [2x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRR1A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3PRR1A1_Jinv: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRR1A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRR1A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:43:36
% EndTime: 2018-12-20 17:43:37
% DurationCPUTime: 0.09s
% Computational Cost: add. (27->21), mult. (54->51), div. (9->3), fcn. (75->14), ass. (0->34)
t33 = koppelP(1,1);
t32 = koppelP(2,1);
t31 = koppelP(3,1);
t30 = koppelP(1,2);
t29 = koppelP(2,2);
t28 = koppelP(3,2);
t27 = xP(3);
t26 = cos(qJ(2,1));
t25 = cos(qJ(2,2));
t24 = cos(qJ(2,3));
t23 = sin(qJ(2,1));
t22 = sin(qJ(2,2));
t21 = sin(qJ(2,3));
t20 = legFrame(1,3);
t19 = legFrame(2,3);
t18 = legFrame(3,3);
t17 = 0.1e1 / t23;
t16 = 0.1e1 / t22;
t15 = 0.1e1 / t21;
t14 = cos(t27);
t13 = sin(t27);
t12 = cos(t20);
t11 = cos(t19);
t10 = cos(t18);
t9 = sin(t20);
t8 = sin(t19);
t7 = sin(t18);
t6 = -t13 * t30 + t14 * t33;
t5 = -t13 * t29 + t14 * t32;
t4 = -t13 * t28 + t14 * t31;
t3 = t13 * t33 + t14 * t30;
t2 = t13 * t32 + t14 * t29;
t1 = t13 * t31 + t14 * t28;
t34 = [(t12 * t26 - t9 * t23) * t17 (t23 * t12 + t9 * t26) * t17 ((t6 * t12 + t9 * t3) * t23 - (t3 * t12 - t9 * t6) * t26) * t17; (t11 * t25 - t8 * t22) * t16 (t22 * t11 + t8 * t25) * t16 ((t5 * t11 + t8 * t2) * t22 - (t2 * t11 - t8 * t5) * t25) * t16; (t10 * t24 - t7 * t21) * t15 (t21 * t10 + t7 * t24) * t15 ((t7 * t1 + t4 * t10) * t21 - (t1 * t10 - t7 * t4) * t24) * t15;];
Jinv  = t34;
