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

function Jinv = P3PRRRR1G1P3A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(2,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G1P3A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G1P3A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G1P3A1_Jinv: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G1P3A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G1P3A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 16:44:25
% EndTime: 2020-08-06 16:44:25
% DurationCPUTime: 0.06s
% Computational Cost: add. (6->6), mult. (18->18), div. (12->6), fcn. (42->18), ass. (0->22)
t19 = cos(qJ(3,3));
t30 = 0.1e1 / t19 / sin(qJ(2,3));
t21 = cos(qJ(3,2));
t29 = 0.1e1 / t21 / sin(qJ(2,2));
t23 = cos(qJ(3,1));
t28 = 0.1e1 / t23 / sin(qJ(2,1));
t27 = t19 * cos(qJ(2,3));
t26 = t21 * cos(qJ(2,2));
t25 = t23 * cos(qJ(2,1));
t18 = sin(qJ(3,1));
t17 = sin(qJ(3,2));
t16 = sin(qJ(3,3));
t15 = legFrame(1,3);
t14 = legFrame(2,3);
t13 = legFrame(3,3);
t6 = cos(t15);
t5 = cos(t14);
t4 = cos(t13);
t3 = sin(t15);
t2 = sin(t14);
t1 = sin(t13);
t7 = [(t3 * t18 + t6 * t25) * t28, (-t6 * t18 + t3 * t25) * t28, 1; (t2 * t17 + t5 * t26) * t29, (-t5 * t17 + t2 * t26) * t29, 1; (t1 * t16 + t4 * t27) * t30, (t1 * t27 - t4 * t16) * t30, 1;];
Jinv  = t7;
