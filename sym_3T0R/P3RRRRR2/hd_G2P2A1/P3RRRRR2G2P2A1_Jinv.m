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

function Jinv = P3RRRRR2G2P2A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(2,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G2P2A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G2P2A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G2P2A1_Jinv: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G2P2A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G2P2A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:39:15
% EndTime: 2020-08-07 03:39:15
% DurationCPUTime: 0.11s
% Computational Cost: add. (15->12), mult. (30->30), div. (24->7), fcn. (66->27), ass. (0->29)
t20 = sin(qJ(2,3));
t25 = cos(qJ(3,3));
t37 = (sin(qJ(1,3)) * cos(qJ(2,3)) + cos(qJ(1,3)) * t20) * t25;
t22 = sin(qJ(2,2));
t26 = cos(qJ(3,2));
t36 = (sin(qJ(1,2)) * cos(qJ(2,2)) + cos(qJ(1,2)) * t22) * t26;
t24 = sin(qJ(2,1));
t27 = cos(qJ(3,1));
t35 = t27 * (sin(qJ(1,1)) * cos(qJ(2,1)) + cos(qJ(1,1)) * t24);
t28 = 0.1e1 / pkin(1);
t34 = 0.1e1 / t20 * t28;
t33 = 0.1e1 / t22 * t28;
t32 = 0.1e1 / t24 * t28;
t31 = 0.1e1 / t25 * t34;
t30 = 0.1e1 / t26 * t33;
t29 = 0.1e1 / t27 * t32;
t23 = sin(qJ(3,1));
t21 = sin(qJ(3,2));
t19 = sin(qJ(3,3));
t18 = legFrame(1,2);
t17 = legFrame(2,2);
t16 = legFrame(3,2);
t9 = cos(t18);
t8 = cos(t17);
t7 = cos(t16);
t6 = sin(t18);
t5 = sin(t17);
t4 = sin(t16);
t1 = [(t6 * t23 + t9 * t35) * t29, (t9 * t23 - t6 * t35) * t29, cos(qJ(1,1) + qJ(2,1)) * t32; (t5 * t21 + t8 * t36) * t30, (t8 * t21 - t5 * t36) * t30, cos(qJ(1,2) + qJ(2,2)) * t33; (t4 * t19 + t7 * t37) * t31, (t7 * t19 - t4 * t37) * t31, cos(qJ(1,3) + qJ(2,3)) * t34;];
Jinv  = t1;
