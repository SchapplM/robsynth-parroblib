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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2020-08-07 10:14
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRRRR7V2G2P2A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(7,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V2G2P2A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V2G2P2A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRRRR7V2G2P2A1_Jinv: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V2G2P2A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V2G2P2A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 10:11:45
% EndTime: 2020-08-07 10:11:45
% DurationCPUTime: 0.05s
% Computational Cost: add. (33->15), mult. (24->18), div. (9->3), fcn. (33->18), ass. (0->10)
t1 = 0.1e1 / (pkin(3) * cos(qJ(2,3) + qJ(3,3)) + cos(qJ(2,3)) * pkin(2) + pkin(1));
t12 = t1 * cos(qJ(1,3));
t2 = 0.1e1 / (pkin(3) * cos(qJ(2,2) + qJ(3,2)) + cos(qJ(2,2)) * pkin(2) + pkin(1));
t11 = t2 * cos(qJ(1,2));
t3 = 0.1e1 / (pkin(3) * cos(qJ(2,1) + qJ(3,1)) + cos(qJ(2,1)) * pkin(2) + pkin(1));
t10 = t3 * cos(qJ(1,1));
t6 = legFrame(1,2);
t5 = legFrame(2,2);
t4 = legFrame(3,2);
t7 = [cos(t6) * t10, -sin(t6) * t10, -sin(qJ(1,1)) * t3; cos(t5) * t11, -sin(t5) * t11, -sin(qJ(1,2)) * t2; cos(t4) * t12, -sin(t4) * t12, -sin(qJ(1,3)) * t1;];
Jinv  = t7;
