% Analytische Jacobi-Matrix für parallelen Roboter
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorpose und aktiven Gelenkkoordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Input:
% xP [4x1]
%   Generalized platform coordinates
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% Jinv [4x4]
%   Analytische Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 17:29
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P4RRRRR2G1P1A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(2,1),zeros(4,3),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4RRRRR2G1P1A1_Jinv: qJ has to be [3x4] (double)');
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4RRRRR2G1P1A1_Jinv: xP has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4RRRRR2G1P1A1_Jinv: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4RRRRR2G1P1A1_Jinv: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4RRRRR2G1P1A1_Jinv: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 17:29:45
% EndTime: 2020-08-07 17:29:46
% DurationCPUTime: 0.13s
% Computational Cost: add. (44->16), mult. (8->24), div. (36->9), fcn. (40->28), ass. (0->15)
t18 = 0.1e1 / pkin(1);
t22 = t18 / sin(qJ(2,4));
t21 = 0.1e1 / sin(qJ(2,3)) * t18;
t20 = 0.1e1 / sin(qJ(2,2)) * t18;
t19 = 0.1e1 / sin(qJ(2,1)) * t18;
t8 = qJ(1,1) + qJ(2,1) + legFrame(1,3);
t7 = qJ(1,2) + qJ(2,2) + legFrame(2,3);
t6 = qJ(1,3) + qJ(2,3) + legFrame(3,3);
t5 = qJ(1,4) + qJ(2,4) + legFrame(4,3);
t17 = xP(4);
t4 = -t17 + t8;
t3 = -t17 + t7;
t2 = -t17 + t6;
t1 = -t17 + t5;
t9 = [cos(t8) * t19, sin(t8) * t19, sin(qJ(3,1)) / cos(qJ(3,1)) * t19, (-koppelP(1,2) * cos(t4) + koppelP(1,1) * sin(t4)) * t19; cos(t7) * t20, sin(t7) * t20, sin(qJ(3,2)) / cos(qJ(3,2)) * t20, (-koppelP(2,2) * cos(t3) + koppelP(2,1) * sin(t3)) * t20; cos(t6) * t21, sin(t6) * t21, sin(qJ(3,3)) / cos(qJ(3,3)) * t21, (-koppelP(3,2) * cos(t2) + koppelP(3,1) * sin(t2)) * t21; cos(t5) * t22, sin(t5) * t22, sin(qJ(3,4)) / cos(qJ(3,4)) * t22, (koppelP(4,1) * sin(t1) - koppelP(4,2) * cos(t1)) * t22;];
Jinv  = t9;
