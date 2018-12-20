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
% Datum: 2018-12-20 17:45
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Jinv = P3PRR1A2_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(4,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PRR1A2_Jinv: qJ has to be [2x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRR1A2_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3PRR1A2_Jinv: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRR1A2_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRR1A2_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:45:01
% EndTime: 2018-12-20 17:45:01
% DurationCPUTime: 0.08s
% Computational Cost: add. (15->15), mult. (24->27), div. (18->4), fcn. (39->11), ass. (0->19)
t22 = 0.1e1 / pkin(2);
t25 = t22 / sin(qJ(2,3));
t24 = 0.1e1 / sin(qJ(2,2)) * t22;
t23 = 0.1e1 / sin(qJ(2,1)) * t22;
t19 = koppelP(3,1);
t16 = koppelP(3,2);
t15 = xP(3);
t14 = legFrame(1,3);
t13 = legFrame(2,3);
t12 = legFrame(3,3);
t8 = cos(t15);
t7 = sin(t15);
t6 = cos(t14);
t5 = cos(t13);
t4 = cos(t12);
t3 = sin(t14);
t2 = sin(t13);
t1 = sin(t12);
t9 = [-t6 * t23, -t3 * t23 ((-t3 * t8 + t6 * t7) * koppelP(1,1) + (t3 * t7 + t6 * t8) * koppelP(1,2)) * t23; -t5 * t24, -t2 * t24 ((-t2 * t8 + t5 * t7) * koppelP(2,1) + (t2 * t7 + t5 * t8) * koppelP(2,2)) * t24; -t4 * t25, -t1 * t25 ((t16 * t8 + t19 * t7) * t4 + (t16 * t7 - t19 * t8) * t1) * t25;];
Jinv  = t9;
