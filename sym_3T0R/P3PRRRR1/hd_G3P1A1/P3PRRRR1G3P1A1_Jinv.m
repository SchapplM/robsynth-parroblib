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
% Datum: 2020-08-06 16:44
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3PRRRR1G3P1A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(2,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G3P1A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G3P1A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G3P1A1_Jinv: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G3P1A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G3P1A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 16:44:44
% EndTime: 2020-08-06 16:44:44
% DurationCPUTime: 0.05s
% Computational Cost: add. (6->6), mult. (12->21), div. (12->6), fcn. (39->18), ass. (0->19)
t18 = cos(qJ(2,1));
t17 = cos(qJ(2,2));
t16 = cos(qJ(2,3));
t15 = sin(qJ(2,1));
t14 = sin(qJ(2,2));
t13 = sin(qJ(2,3));
t12 = legFrame(1,2);
t11 = legFrame(2,2);
t10 = legFrame(3,2);
t9 = 0.1e1 / t15;
t8 = 0.1e1 / t14;
t7 = 0.1e1 / t13;
t6 = cos(t12);
t5 = cos(t11);
t4 = cos(t10);
t3 = sin(t12);
t2 = sin(t11);
t1 = sin(t10);
t19 = [(t15 * t3 + t6 * t18) * t9, (t6 * t15 - t3 * t18) * t9, sin(qJ(3,1)) / cos(qJ(3,1)) * t9; (t14 * t2 + t5 * t17) * t8, (t5 * t14 - t2 * t17) * t8, sin(qJ(3,2)) / cos(qJ(3,2)) * t8; (t13 * t1 + t4 * t16) * t7, (-t1 * t16 + t4 * t13) * t7, sin(qJ(3,3)) / cos(qJ(3,3)) * t7;];
Jinv  = t19;
