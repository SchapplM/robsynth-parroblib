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

function Jinv = P3RRRRR2G1P1A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(2,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G1P1A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G1P1A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G1P1A1_Jinv: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G1P1A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G1P1A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:38:58
% EndTime: 2020-08-07 03:38:58
% DurationCPUTime: 0.04s
% Computational Cost: add. (12->6), mult. (0->9), div. (21->7), fcn. (21->15), ass. (0->8)
t7 = 0.1e1 / pkin(1);
t10 = 0.1e1 / sin(qJ(2,3)) * t7;
t9 = 0.1e1 / sin(qJ(2,2)) * t7;
t8 = 0.1e1 / sin(qJ(2,1)) * t7;
t3 = qJ(1,1) + qJ(2,1) + legFrame(1,3);
t2 = qJ(1,2) + qJ(2,2) + legFrame(2,3);
t1 = qJ(1,3) + qJ(2,3) + legFrame(3,3);
t4 = [cos(t3) * t8, sin(t3) * t8, sin(qJ(3,1)) / cos(qJ(3,1)) * t8; cos(t2) * t9, sin(t2) * t9, sin(qJ(3,2)) / cos(qJ(3,2)) * t9; cos(t1) * t10, sin(t1) * t10, sin(qJ(3,3)) / cos(qJ(3,3)) * t10;];
Jinv  = t4;
