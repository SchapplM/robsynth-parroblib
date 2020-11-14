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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d3]';
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
% Datum: 2019-05-03 14:40
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3PPR1G1P1A2_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PPR1G1P1A2_Jinv: qJ has to be [2x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PPR1G1P1A2_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PPR1G1P1A2_Jinv: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PPR1G1P1A2_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PPR1G1P1A2_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:40:03
% EndTime: 2019-05-03 14:40:03
% DurationCPUTime: 0.04s
% Computational Cost: add. (9->9), mult. (18->18), div. (0->0), fcn. (24->8), ass. (0->19)
t18 = koppelP(1,1);
t17 = koppelP(2,1);
t16 = koppelP(3,1);
t15 = koppelP(1,2);
t14 = koppelP(2,2);
t13 = koppelP(3,2);
t12 = xP(3);
t11 = legFrame(1,3);
t10 = legFrame(2,3);
t9 = legFrame(3,3);
t8 = cos(t12);
t7 = sin(t12);
t6 = cos(t11);
t5 = cos(t10);
t4 = cos(t9);
t3 = sin(t11);
t2 = sin(t10);
t1 = sin(t9);
t19 = [t6, t3, (-t15 * t3 - t18 * t6) * t7 + (-t15 * t6 + t18 * t3) * t8; t5, t2, (-t14 * t2 - t17 * t5) * t7 + (-t14 * t5 + t17 * t2) * t8; t4, t1, (-t1 * t13 - t16 * t4) * t7 + (t1 * t16 - t13 * t4) * t8;];
Jinv  = t19;
