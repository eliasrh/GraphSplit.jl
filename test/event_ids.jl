@testset "exact event-ID input range" begin
    mktempdir() do directory
        catalog=joinpath(directory,"catalog.txt")
        maxid="9007199254740991"
        write(catalog,"2026 1 1 0 0 0 64 -20 3 1 $maxid\n")
        @test only(GraphSplit.read_catalog(catalog).event_id)==parse(Int64,maxid)
        for id in ("9007199254740992","9007199254740993","999999999999999999999999999", "3.5", "0", "-1")
            write(catalog,"2026 1 1 0 0 0 64 -20 3 1 $id\n")
            @test_throws ErrorException GraphSplit.read_catalog(catalog)
        end
        theta=joinpath(directory,"theta");mkpath(theta)
        path=joinpath(theta,"theta_STA_P.txt")
        write(path,"$maxid 0.1 7\n7 0 7\n")
        group=only(GraphSplit.load_theta_folder(theta))
        @test haskey(group.entries,parse(Int64,maxid))
        for text in ("9007199254740993 0.1 7\n", "7 0.1 9007199254740993\n", "7.5 0.1 7\n", "7 0.1 7.5\n")
            write(path,text)
            @test_throws ErrorException GraphSplit.load_theta_folder(theta)
        end
        std=joinpath(directory,"std.txt")
        write(std,"9007199254740993 0.002 7 10\n")
        @test_throws ErrorException GraphSplit.read_std_entries(std)
        shift=joinpath(directory,"shift.txt")
        write(shift,"0 0 0 0 9007199254740993\n")
        @test_throws ErrorException GraphSplit.read_catalog_shift(shift,Int64[7])
        # Shuffling catalog rows preserves the ID-to-row correspondence.
        write(catalog,"2026 1 1 0 0 0 64 -20 3 1 81\n2026 1 1 0 0 0 64.1 -20 3 1 7\n")
        mapping=GraphSplit.event_row_map(GraphSplit.read_catalog(catalog))
        @test mapping[81]==1
        @test mapping[7]==2
    end
end
