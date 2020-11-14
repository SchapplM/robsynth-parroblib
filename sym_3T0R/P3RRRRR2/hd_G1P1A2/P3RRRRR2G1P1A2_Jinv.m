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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2020-08-07 03:39
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRRRR2G1P1A2_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(2,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G1P1A2_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G1P1A2_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G1P1A2_Jinv: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G1P1A2_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G1P1A2_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:39:06
% EndTime: 2020-08-07 03:39:06
% DurationCPUTime: 0.14s
% Computational Cost: add. (78->39), mult. (36->33), div. (30->11), fcn. (45->36), ass. (0->21)
t28 = -2 * pkin(1);
t27 = 1 / pkin(2) / pkin(1);
t12 = qJ(1,1) + legFrame(1,3);
t11 = qJ(1,2) + legFrame(2,3);
t10 = qJ(1,3) + legFrame(3,3);
t26 = 0.1e1 / (sin(qJ(2,3) + qJ(3,3)) + sin(qJ(2,3) - qJ(3,3))) * t27;
t25 = 0.1e1 / (sin(qJ(2,2) + qJ(3,2)) + sin(qJ(2,2) - qJ(3,2))) * t27;
t24 = 0.1e1 / (sin(qJ(2,1) + qJ(3,1)) + sin(qJ(2,1) - qJ(3,1))) * t27;
t23 = qJ(2,1) + t12;
t22 = qJ(2,2) + t11;
t21 = qJ(2,3) + t10;
t18 = cos(qJ(3,1));
t17 = cos(qJ(3,2));
t16 = cos(qJ(3,3));
t9 = -qJ(3,1) + t23;
t8 = qJ(3,1) + t23;
t7 = -qJ(3,2) + t22;
t6 = qJ(3,2) + t22;
t5 = -qJ(3,3) + t21;
t4 = qJ(3,3) + t21;
t1 = [(cos(t12) * t28 + (-cos(t9) - cos(t8)) * pkin(2)) * t24, (sin(t12) * t28 + (-sin(t9) - sin(t8)) * pkin(2)) * t24, -(cos(qJ(2,1)) * pkin(1) + t18 * pkin(2)) * sin(qJ(3,1)) / t18 ^ 2 / sin(qJ(2,1)) * t27; (cos(t11) * t28 + (-cos(t7) - cos(t6)) * pkin(2)) * t25, (sin(t11) * t28 + (-sin(t7) - sin(t6)) * pkin(2)) * t25, -(cos(qJ(2,2)) * pkin(1) + t17 * pkin(2)) * sin(qJ(3,2)) / t17 ^ 2 / sin(qJ(2,2)) * t27; (cos(t10) * t28 + (-cos(t5) - cos(t4)) * pkin(2)) * t26, (sin(t10) * t28 + (-sin(t5) - sin(t4)) * pkin(2)) * t26, -(cos(qJ(2,3)) * pkin(1) + t16 * pkin(2)) * sin(qJ(3,3)) / t16 ^ 2 / sin(qJ(2,3)) * t27;];
Jinv  = t1;
